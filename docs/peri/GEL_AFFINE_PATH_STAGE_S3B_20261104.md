# S3b material-point affine-F path and consistent local tangent

Research candidate only, not production PD integration or a complete S3 stage.
For F(s)=F0+s*dF, deviatoric Green strain is E0+sG1+s²G2, with
G1=dev sym(F0^T*dF), G2=dev(dF^T*dF)/2. Maxwell extension X obeys
X'+hX=G1+2sG2, h=dt/tau. The analytic solution is

    X(s)=exp(-hs)X0 + (1-exp(-hs))/h*G1
         +2*(s/h-(1-exp(-hs))/h²)*G2.

End history is analytic, not a linear interpolation of endpoint Green strain.
Physical loss is2Gv*h*int||X||² ds. Native uses analytic3x3 Gram moments with
small-h series; negative increments fail without clipping. Averaged PK1 uses
fixed32-point Gauss-Legendre on this same actual path, with the unchanged
neo-Hookean/volumetric and Green-Maxwell constitutive energy. No stress/work
rescaling enforces the ledger. End storage + physical loss must independently
match int PK1:dF to normalized1e-10 before history commit. Cubic determinant
stationary points guard the whole interval against J<.05, not only endpoints.
Zero dt permits only unchanged F and returns the complete accepted history.

The9x9 endpoint tangent analytically differentiates average PK1, including the
path/history terms. It is recorded as column-major operator blocks, each block
itself row-major3x3. A causal viscous path can yield a nonsymmetric consistent
tangent. Symmetry is measured, never manufactured; use of an SPD-only global
solver therefore requires a separate admissibility/algorithm derivation.
Zero branch has the conservative symmetric averaged tangent. This is an
operator qualification issue, not a failed material arithmetic gate to waive.

Independent reference uses65-digit Decimal polynomial-plus-exponential
particular solutions and dimensional closed antiderivatives for loss, rather
than native Gram/series. Stress uses48-point high-precision quadrature checked
against polynomial moments; the9x9 tangent is independently differentiated
with65-digit central differences at epsilon1e-18. All normalized fields have
gate1e-10. Superposed constant rotations, relaxation, ramp/hold/release, zero
branch, volumetric and isochoric shear/mixed paths, prehistory and no advance
are required. Negative controls must reject wrong linear-Green history and
symmetrization of the true local tangent.

Freeze source/ELF/libs/plan/inputs/reader before native dispatch; each native5s,
256MiB address space/RSS,64 path segments,32MiB raw. F host/output reserves stay
120/150GiB. Every native PID/parent/argv/startticks/boot and exact return code/
raw/size/SHA is retained. Unexpected execution/integrity/resource/accuracy stops
later cases without fallback/reset or redispatching an old state. A separate
single immutable review is required; runtime checks are not release.

Global PD work-conjugate gradient/transpose, nonlinear constrained solver and
its admissible nonsymmetric/symmetric-part properties, time/space/horizon/
large-strain qualification, full production history restart and experimental
calibration remain required after this material-point control.
