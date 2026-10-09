# FLOW36_EM: GPU optimization log

Changes made to FLOW36_EM (polymer version of FLOW36) to fix GPU performance and multi-node
scaling on Leonardo (4x A100 64 GB per node, NVHPC + CUDA-aware MPI). Baseline: commit `054bf76`.
This file is the reference for porting the same changes to the official FLOW36 repository
(GPU branch).

## Results

256^3, single phase (Leonardo, 4 GPUs per node), time per step (steps 3-9):

| Version                                                    | 1 node (4 GPUs) | 2 nodes (8 GPUs)   |
|------------------------------------------------------------|-----------------|--------------------|
| Original code                                              | 90 ms           | slower than 1 node |
| Transposes: one exchange, persistent buffers               | ~85 ms          | ~80 ms             |
| + `calculate_var` boundary sums in one kernel              | ~32 ms          |                    |
| + fused `helmholtz` solver                                 | ~29-30 ms       |                    |
| + device (cudaMalloc) transpose buffers, GPUDirect RDMA    | ~26 ms          | ~29 ms             |
| + persistent FFT/DCT/transform work arrays, no copies      | **24 ms**       | **25-26 ms**       |

512^3, phase field + polymers (Saramito, `fenep_flag=3`, `pol_res_flag=0`):

| Version                                                    | 4 nodes (16 GPUs) | 8 nodes (32 GPUs) |
|------------------------------------------------------------|-------------------|-------------------|
| Before the polymer/phase-field fixes                       | ~60 s             |                   |
| After (batched polymer Helmholtz, phase-field loops on GPU)| **1.3 s**         | **0.77 s** (84% efficiency) |

GPU memory, 512^3 phase field + polymers on 32 GPUs: 5.6 GB per GPU (printed after step 1).
Correctness: results checked against the original code (diagnostics agree); the rewritten
Helmholtz solvers are bitwise identical to the original ones on CPU.

## What was wrong (root causes)

1. **Transposes exchanged with one partner at a time** (P-1 sequential rounds of
   pack/isend/recv/wait/unpack, plus a useless `bufs=0*bufs` kernel per round).
2. **Allocation at every call with `-gpu=managed`**: every `allocate` of a managed array is
   page-faulted onto the GPU at first touch, every `deallocate` synchronizes the device. The
   transforms, FFT/DCT routines and transposes allocated/freed >100 arrays per time step.
3. **Managed memory in MPI buffers**: the network card cannot access managed memory (no
   GPUDirect RDMA), so inter-node messages were staged through the host (~0.3 GB/s).
4. **Loops running on the host** (implicit, no warning at run time):
   - array-syntax loops over `k` inside `!$acc kernels` launching 4 tiny kernels per `k`
     (`calculate_var`, ~1000 launches per wall per call);
   - loop nests over full fields with no OpenACC directive (105 in `phi_non_linear`,
     `sterm_ch`, `calculate_var`): whole fields migrate GPU<->host at every loop;
   - `calculate_pol`: host loop over all (i,j) columns calling the 1D solver `helmholtz_rred`
     for every conformation component (~1e5 calls, ~1e6 kernels per step at 512^3).
5. **Helmholtz solver**: assembly of 5 full 3D coefficient arrays in 4 kernels + solve + copy,
   for coefficients that are simple formulas of `beta2(i,j)` and `k`.

## Code changes (`source_code/`)

### Transposes: `xy2xz.f90`, `xz2xy.f90`, `xz2yz.f90`, `yz2xz.f90`, `module.f90`, `main.f90`
Commits `9761eaf`, `9ec53db`, `32263f3`, `23a0d57`.
- Each transpose (and its `_fg` variant) is a thin wrapper calling one core routine `*_a2a`
  with the global sizes (`nx/ny/nz` or `npsix/npsiy/npsiz`): one pack kernel for all the
  blocks, all `mpi_irecv`/`mpi_isend` posted at once + one `mpi_waitall`, one unpack kernel.
  The own block (j=me) is not sent: the unpack kernel takes it from the send buffer.
  Not `mpi_alltoall`: with managed memory the collective copied the own block on the host.
- Blocks keep the original padded size `ngx*ngz*ngy*2`; block j = partner with coordinate j.
- `main.f90`: sub-communicators `cart_comm_dir(0:1)` (`mpi_cart_sub` of `cart_comm`), with a
  startup check that the sub-communicator rank equals the Cartesian coordinate; freed at the end.
- `module.f90`, module `a2a_buffers`: `a2a_send`/`a2a_recv` allocated once (`a2a_reserve`),
  **CUDA Fortran `device` arrays** under `#ifdef _CUDA` (ordinary arrays otherwise), indexed 1D
  in the kernels (block j at offset `j*numel`), passed directly to the MPI calls.

### Transforms and FFT/DCT: `spectral_to_phys.f90`, `phys_to_spectral.f90`, `fftx_*`, `ffty_*`, `dctz_*`
Commit `fe6e5ca` (normal grid only; `_fg` variants unchanged).
- FFT/DCT work arrays (`wt`, `wot`, `a`, `b`, `ac`, `bc`) are `allocatable, save`, reallocated
  only if the shape changes; no `deallocate` at every call.
- `spectral_to_phys`/`phys_to_spectral`: one persistent array per stage (`s1`,`s2`,`s3`,
  allocated with `stage_alloc` from `a2a_buffers`) and a pointer to the current stage; the
  transposes write directly into the next stage, removing the `wa=u` copies. The `u=uc` copy
  before `dctz_bwd` is kept (it modifies its input).

### Boundary sums: `calculate_var.f90`
Commit `7fe84e3`. The four blocks computing `a11,a12,fr1,fc1` / `a21,a22,fr2,fc2` (no-slip and
free-slip at both walls): one `parallel loop collapse(2)` over (i,j) with the sum over k in
registers, same terms in the same order.

### Helmholtz: `helmholtz.f90`, `calculate_pol.f90`
Commits `871ec6e`, `9ce9afa`.
- `helmholtz`: assembly + `gauss_solver` + `f=h` fused in one kernel per (i,j) column:
  `a,c` and initial `b` computed on the fly, `d,e` as 1D arrays, values modified during the
  elimination kept in parity slots (used two steps later); only the final RHS and diagonal are
  stored (persistent). Bitwise identical on CPU (Dirichlet/Neumann/mixed, nz=257 and 64).
- `helmholtz_cols`: same kernel with boundary values per column, `r1(i,j,c)`, `r2(i,j,c)`.
- `calculate_pol`: 6 calls of `helmholtz_cols` (one per conformation component) instead of the
  host loop over columns calling `helmholtz_rred`. Bitwise identical on CPU. `gauss_solver`
  kept for `helmholtz_pol`; `helmholtz_fg`, `helmholtz_red`, `helmholtz_pol` unchanged.

### Host loops to GPU: `phi_non_linear.f90`, `sterm_ch.f90`, `calculate_var.f90`
Commit `9ce9afa`. 105 element-wise loop nests got `!$acc parallel loop collapse(n)` (checked:
no calls/IO, no reads of the written array at other indices).

### Build fixes needed with NVHPC 24.3/25.11 (pre-existing issues)
- `helmholtz.f90` (`helmholtz_pol`): `!$acc loop` instead of `parallel loop` nested in an open
  `kernels` region (NVFORTRAN-S-0155). Commit `8838b48`.
- `assemble.f90` (`assemble_eigen`): loop calling `Diagonalization` (LAPACK `DGEEV`, allocate)
  runs on the host; LAPACK/BLAS linked. Commit `2d05cd2`.
- `initialize_pol.f90`: `assemble_eigen` on restart only with `solvpolflag==1` (its arrays are
  allocated only then). Commit `527af4e`.
- `courant_check.f90`: NaN checks on host copies (`chk_phi`, `chk_psi`) of `phic(1,1,1,1)`,
  `psic_fg(1,1,1,1)`: with `-cuda` + managed memory `ieee_is_nan` on a managed element resolved
  to the device routine (NVFORTRAN-S-0155). Commit `8e75842`.

### Diagnostics: `main.f90`
- NVTX range per time step (`step <n>`) under `#ifdef USE_NVTX` (commit `6e247b0`).
- GPU memory in use after the first step, max over ranks, under `#ifdef _CUDA` (`34f37c1`).

## Build and run setup (Leonardo, `machine="22"`)

- `compile.sh` (machine 22): `nvhpc/25.11` + HPC-X 2.25.1 bundled with NVHPC
  (`comm_libs/12.9/hpcx/hpcx-2.25.1/hpcx-init.sh`, `hpcx_load`); copies `binder_leo.sh` to
  `set_run`. Commits `4cab023`, `d01f4ff`.
- `Leonardo/makefile_gpu`: `FC=mpif90`; `-fast -acc -cuda -gpu=managed,cuda12.9`; no
  `-Minfo=accel`, no `-lnvToolsExt` (removed in CUDA 12.9); `LIBS=-cudalib=cufft -llapack -lblas`;
  `NVTX=1` adds `-DUSE_NVTX -cudalib=nvtx`.
- `Leonardo/go_gpu.sh`: `--ntasks=NUMTASKS` (job sized by NYCPU*NZCPU); same modules as the
  build; `cd` into `set_run` (submit with `sbatch set_run/go.sh` from the main folder);
  `mpirun --map-by ppr:4:node --mca pml ucx ./binder_leo.sh ./sc_compiled/flow36`;
  `PROFILE=1` profiles every rank with nsys (`report_<rank>.nsys-rep`).
- `Leonardo/binder_leo.sh`: `UCX_NET_DEVICES=mlx5_<local rank>:1` (from MHIT36).
- Process grid: `NYCPU=4` so that each node holds one y-group (xz2yz/yz2xz stay on NVLink).

## Porting to the official FLOW36 GPU branch (suggested order)

Test after each step (single phase 256^3, 1 and 2 nodes, compare diagnostics with the original).

1. Build setup: makefile (`mpif90`, `-cuda`, `cuda12.x`, no `-lnvToolsExt`), compile.sh
   machine block, job script, binder. Fix the pre-existing build errors (section above).
2. Transposes + `a2a_buffers` + `cart_comm_dir` (biggest multi-node effect; device buffers).
3. Transforms/FFT/DCT persistent arrays.
4. `calculate_var` boundary sums; fused `helmholtz`.
5. If the official version has the polymer/phase-field code: `helmholtz_cols` + `calculate_pol`,
   GPU directives on the host loops (search loop nests over `fpy/spy` without `!$acc`).
6. NVTX and memory report.

Things to check in the official version: whether it has the same host loops (compile once with
`-Minfo=accel` and grep `Sequential loop scheduled on host`, `serial kernel`), and the same
`allocate`/`deallocate` pattern in the transforms.

## Open issues (not changed, to be decided)

- `phys_to_spectral(u,uout)` has no `aliasing` argument but is called with 3 arguments:
  `aliasing` is an uninitialized local, so dealiasing of the forward transform is undefined.
  Fix: add the argument (as in `phys_to_spectral_fg`); this can change results.
- `calculate_pol`: the imaginary part of `Cxz` uses the boundary value of the real part
  (`bc_cxz(i,nz,j,1)`), because `s_bc(1)=bc_cxz(i,nz,j,2)` was commented out in the original.
  Reproduced as is.
- `print_start.f90`: Bingham number printed with `f8.5` (`********` for `Bi=200`).
- Not yet optimized: `_fg` (dual grid) transforms and solvers, `helmholtz_red`, asynchronous
  OpenACC regions (one `cuStreamSynchronize` per region), batching several fields per transform.
