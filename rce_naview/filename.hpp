#pragma once

#include <cctype>
#include <string>

// Stem used for the SVG file. Characters outside [A-Za-z0-9._] become '_',
// then the stem is cut at 80 characters. Empty, dot, dot-dot, and Windows
// device names are replaced with record_<line_no>.
inline bool is_windows_device_name(const std::string& stem) {
    std::string u;
    u.reserve(stem.size());
    for (unsigned char c : stem) {
        u.push_back(static_cast<char>(std::toupper(c)));
    }
    if (u == "CON" || u == "PRN" || u == "AUX" || u == "NUL") {
        return true;
    }
    if (u.size() == 4 && (u.compare(0, 3, "COM") == 0 || u.compare(0, 3, "LPT") == 0)) {
        char d = u[3];
        return d >= '1' && d <= '9';
    }
    return false;
}

inline std::string stem_from_label(const std::string& label, int line_no) {
    std::string stem;
    stem.reserve(label.size());
    for (unsigned char c : label) {
        bool keep = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                    (c >= '0' && c <= '9') || c == '.' || c == '_';
        stem.push_back(keep ? static_cast<char>(c) : '_');
    }
    if (stem.size() > 80) {
        stem.resize(80);
    }
    if (stem.empty() || stem == "." || stem == ".." || is_windows_device_name(stem)) {
        stem = "record_" + std::to_string(line_no);
    }
    return stem;
}
