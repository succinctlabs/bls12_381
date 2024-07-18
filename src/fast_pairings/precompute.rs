use crate::{fp::Fp, fp2::Fp2, G1Affine, Scalar};

pub(crate) fn line_evaluation(p: &G1Affine, alpha: Fp2, bias: Fp2) -> (Fp2, <Fp as Mul<Fp2>>::Output, Fp2) {
    (
        -bias,
        -Scalar::from(p.x) * alpha,
        Fp2 {
            Fp::zero(),
            p.y
        }
    )
}
