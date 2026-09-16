# Beyond univariate analysis: predicting species richness with model-based ordination

Lunchtime workshop at the [Joint NSO–GfÖ Conference 2026](https://nordicsocietyoikos.glueup.com/event/joint-nso-gf%c3%b6-conference-2026-135623/workshops.html#workshops), 13–17 September 2026, Odeon, Odense.

**When:** Thursday 17 September 2026, 12:30–14:00

**Organiser:** [Bert van der Veen](https://bertvanderveen.github.io/) (bert.van-der-veen@uni-bayreuth.de)

## About

Species richness is one of the most commonly analysed metrics in community ecology, calculated as the sum of species presences at sites. Typically, it is analysed using a univariate regression, which violates various assumptions, including the mean-variance relationship, independence, and heterogeneity of species responses, with potential for erroneous conclusions. In this session, I will demonstrate how species richness can be analysed using multispecies models in the `gllvm` R-package.

## Materials

- [`presentation_short.pdf`](presentation_short.pdf): introductory slides (source in [`presentation/`](presentation/))
- [`practical.Rmd`](practical.Rmd) / [`practical.pdf`](practical.pdf): the practical, analysing species richness of the Skabbholmen vegetation data with `gllvm`

## Requirements

The most recent version of R and the `gllvm` package:

```r
install.packages("gllvm")
```

## Additional resources

- [`gllvm` vignettes](https://jenniniku.github.io/gllvm/): package documentation and worked examples
- [GLLVM workshop](https://github.com/BertvanderVeen/GLLVM-workshop): material from the Physalia course on generalized linear latent variable models
- Rugstad, A., O'Hara, B., van der Veen, B. & Mehlhoop, A. C. (2026). One Toolbox, Many Tools: A Practitioner's Guide to Model-based Ordination for Community Ecology. *EcoEvoRxiv* (preprint). <https://doi.org/10.32942/X2KM2V>
