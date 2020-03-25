abstract type AbstractFPGroup <: Group end

@doc doc"""
    ::GSymbol
> Represents a syllable.
> Abstract type which all group symbols of AbstractFPGroups should subtype. Each
> concrete subtype should implement fields:
> * `id` which is the `Symbol` representation/identification of a symbol
> * `pow` which is the (multiplicative) exponent of a symbol.
"""
abstract type GSymbol end
