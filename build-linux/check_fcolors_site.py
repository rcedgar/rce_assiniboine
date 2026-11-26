from pathlib import Path

root = Path("/mnt/c/src/rce_assiniboine")
fcolors = root / "big_control_website_fcolors"
default = root / "big_control_website"
index = (fcolors / "index.html").read_text(encoding="utf-8")
print("fcolors thumb", "big_control_thumbnail_pngs_w20/" in index)
print("fcolors rowbreak", index.count("rowbreak"))
print("fcolors class vle", index.count('class="vle"'))
print("fcolors inline background", index.count('style="background'))
zet = (fcolors / "VLE_zetavirus_index.html").read_text(encoding="utf-8")
ane = (fcolors / "DECOY_anellovirus_index.html").read_text(encoding="utf-8")
print("zet class", "thumbs vle" in zet, 'class="thumbs"' in zet)
print("ane class", "thumbs decoy" in ane, 'class="thumbs"' in ane)
seq = (
    fcolors
    / "VLE_ambivirus_SerratusID90_SRR8549884_NODE_80_length_3141_cov_6.603651.html"
).read_text(encoding="utf-8")
print("seq w20", "big_control_svgs_100_per_family_w20/" in seq)
print("seq rainbow dir", "big_control_svgs_100_per_family/" in seq)
svg = (
    root
    / "big_control_svgs_100_per_family_w20"
    / "VLE_ambivirus_SerratusID90_SRR8549884_NODE_80_length_3141_cov_6.603651.svg"
).read_text(encoding="utf-8")
print("svg hsl", svg.count("hsl("))
print("svg ffa360", "#ffa360" in svg)
print("svg 8c8c8c", "#8c8c8c" in svg)
dindex = (default / "index.html").read_text(encoding="utf-8")
print("default vle", dindex.count('class="vle"'))
print("default black png", "big_control_thumbnail_pngs/" in dindex)
print("default w20", "w20" in dindex)
