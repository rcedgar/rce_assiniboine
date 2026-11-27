#pragma once

#include <string>

// One circle scored by the published one-term VLE formula
// S = (pct_bulge - 12) / 5.21.
struct VlenessScore {
    int length;
    double pct_bulge;
    double score;
};

// Score a circular dot-bracket of '(' ')' and '.'.
// On success fills score and returns true.
// On a malformed structure returns false and, when error is not null,
// sets a short reason. Does not throw.
bool vleness_score(const std::string& dot_bracket, VlenessScore* score, std::string* error);
