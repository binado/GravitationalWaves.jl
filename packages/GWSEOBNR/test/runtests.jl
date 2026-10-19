using Test, GWSEOBNR
@test approximant(SEOBNRv5HM()) == "SEOBNRv5HM"
@test !is_differentiable(SEOBNRv5HM())
