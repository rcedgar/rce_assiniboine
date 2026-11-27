# Features in the VLE score

S is a weighted sum of standardized features. Each feature has its training median subtracted and is then divided by its training interquartile range. Counts and lengths are replaced by their logarithm before that step.

A minus sign on a term means a lower value votes for a viral-like element. A term without a minus means a higher value votes for one.

Call the circle a viral-like element when S is at least 0.710173. That cutoff is the point with the highest recall whose precision is still at least 0.5.

Formula

S = 1*z(pct_bulge)

## pct_bulge

Weight 1. Sign +1. Higher values vote for a viral-like element.

Unpaired nucleotides inside composite stems, as a percentage of the circle. This is the share of the circle spent in bulges and interior loops rather than in pairs, hairpin loops, or branch points.
