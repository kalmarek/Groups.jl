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

abstract type GWord{T<:GSymbol} <: GroupElem end

@doc doc"""
    W::GroupWord{T} <: GWord{T<:GSymbol} <:GroupElem
> Basic representation of element of a finitely presented group. `W.symbols`
> fieldname contains particular group symbols which multiplied constitute a
> group element, i.e. a word in generators.
> As reduction (inside group) of such word may be time consuming we provide
> `savedhash` and `modified` fields as well:
> hash (used e.g. in the `unique` function) is calculated by reducing the word,
> setting `modified` flag to `false` and computing the hash which is stored in
> `savedhash` field.
> whenever word `W` is changed `W.modified` is set to `false`;
> Future comparisons don't perform reduction (and use `savedhash`) as long as
> `modified` flag remains `false`.

"""

mutable struct GroupWord{T} <: GWord{T}
    symbols::Vector{T}
    modified::Bool
    savedhash::UInt
    parent::Group

    function GroupWord{T}(symbols::Vector{<:GSymbol}) where T
       return new{T}(symbols, true, zero(UInt))
    end
end
