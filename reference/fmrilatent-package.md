# fmrilatent: Latent Space Representations of fMRI Data

\`fmrilatent\` owns latent representations for neuroimaging data,
including explicit factorized objects, shared basis assets,
operator-backed latent objects, coefficient recovery, and
coefficient-to-map projection.

## Details

The package deliberately stops at the latent-method layer. It does not
own first-level GLM fitting, contrasts, temporal autocorrelation
modeling, or statistical inference. Downstream modeling packages should
consume \`coef_time(x, "analysis")\` plus the decoder and projection
generics exposed here.

## See also

Useful links:

- <https://github.com/bbuchsbaum/fmrilatent>

- Report bugs at <https://github.com/bbuchsbaum/fmrilatent/issues>

## Author

**Maintainer**: Bradley Buchsbaum <brad.buchsbaum@gmail.com>
([ORCID](https://orcid.org/0000-0002-1108-4866))
