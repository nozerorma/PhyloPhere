# Palettes

# Comment, uncomment and add families as needed depending on the scope of the analysis.

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Load paletteer for colorblind-friendly palettes
suppressPackageStartupMessages({
  library(paletteer)
})

# Tableau 10 colorblind-friendly palette as base
tableau_10_cb <- paletteer_d("ggthemes::Tableau_10")

# Helper function to darken colors
darken_color <- function(hex, factor = 0.6) {
  rgb_vals <- col2rgb(hex) / 255
  darkened <- rgb_vals * factor
  rgb(darkened[1], darkened[2], darkened[3])
}

## Colorblind-friendly primate color palette
## Using Tableau 10 and Paul Tol palettes for maximum accessibility
## Removed Callitrichidae as NCBI has collapsed it into Cebidae, but can be added back if needed with a distinct color
primates_palette <- c(
  "Aotidae"           = "#4E79A7",  # Blue (Tableau)
  # "Callitrichidae"    = "#BAB0AC",  # Gray (Tableau)
  "Cebidae"           = "#E15759",  # Vermillion (Tableau)
  "Atelidae"          = "#76B7B2",  # Bluish_green (Tableau)
  "Cercopithecidae"   = "#EDC948",  # Yellow (Tableau)
  "Hominidae"         = "#B07AA1",  # Reddish_purple (Tableau)
  "Hylobatidae"       = "#FF9DA7",  # Pink (Tableau)
  "Lemuridae"         = "#F28E2B",  # Orange (Tableau)
  "Indriidae"         = "#9C755F",  # Brown (Tableau)
  "Cheirogaleidae"    = "#59A14F",  # Green (Tableau)
  "Lorisidae"         = "#77B4C8",  # Light blue
  "Galagidae"         = "#C9A961",  # Tan
  "Pitheciidae"       = "#A0796A",  # Earth tone
  "Daubentoniidae"    = "#D37295",  # Mauve
  "Tarsiidae"         = "#8E8DBE"   # Lavender
)

## Darkened primate color palette (60% darker for contrast)
dark_primates_palette <- c(
  "Aotidae"         = darken_color("#4E79A7", 0.6),
  # "Callitrichidae"  = darken_color("#BAB0AC", 0.6),
  "Cebidae"         = darken_color("#E15759", 0.6),
  "Atelidae"        = darken_color("#76B7B2", 0.6),
  "Cercopithecidae" = darken_color("#EDC948", 0.6),
  "Hominidae"       = darken_color("#B07AA1", 0.6),
  "Hylobatidae"     = darken_color("#FF9DA7", 0.6),
  "Lemuridae"       = darken_color("#F28E2B", 0.6),
  "Indriidae"       = darken_color("#9C755F", 0.6),
  "Cheirogaleidae"  = darken_color("#59A14F", 0.6),
  "Lorisidae"       = darken_color("#77B4C8", 0.6),
  "Galagidae"       = darken_color("#C9A961", 0.6),
  "Pitheciidae"     = darken_color("#A0796A", 0.6),
  "Daubentoniidae"  = darken_color("#D37295", 0.6),
  "Tarsiidae"       = darken_color("#8E8DBE", 0.6)
)

# Colorblind-friendly mammalian order color palette
mammals_palette <- c(
  "Carnivora"        = "#4E79A7",  # Blue
  "Perissodactyla"   = "#E15759",  # Vermillion
  "Hyracoidea"       = "#76B7B2",  # Bluish_green
  "Afrosoricida"     = "#59A14F",  # Green
  "Artiodactyla"     = "#EDC948",  # Yellow
  "Rodentia"         = "#F28E2B",  # Orange
  "Chiroptera"       = "#B07AA1",  # Reddish_purple
  "Diprotodontia"    = "#9C755F",  # Brown
  "Scandentia"       = "#BAB0AC",  # Gray
  "Proboscidea"      = "#77B4C8",  # Light blue
  "Sirenia"          = "#C9A961",  # Tan
  "Tubulidentata"    = "#8E8DBE",  # Lavender
  "Macroscelidea"    = "#D37295",  # Mauve
  "Eulipotyphla"     = "#A0796A",  # Earth
  "Lagomorpha"       = "#FF9DA7",  # Pink
  "Pholidota"        = "#7C7C7C",  # Medium gray
  "Dermoptera"       = "#5F5F5F"   # Dark gray
)

## Compatibility aliases
primate_family_colors <- primates_palette
dark_family_palette <- dark_primates_palette
mammal_order_colors <- mammals_palette

## Colorblind-friendly fallback palette using Tableau 20
fallback_palette <- c(
  "Category1"  = "#4E79A7",  # Blue
  "Category2"  = "#F28E2B",  # Orange
  "Category3"  = "#E15759",  # Vermillion
  "Category4"  = "#76B7B2",  # Bluish_green
  "Category5"  = "#59A14F",  # Green
  "Category6"  = "#EDC948",  # Yellow
  "Category7"  = "#B07AA1",  # Reddish_purple
  "Category8"  = "#FF9DA7",  # Pink
  "Category9"  = "#9C755F",  # Brown
  "Category10" = "#BAB0AC",  # Gray
  "Category11" = "#8CD17D",  # Light green
  "Category12" = "#B6992D",  # Olive
  "Category13" = "#499894",  # Teal
  "Category14" = "#86BCB6",  # Sage
  "Category15" = "#F1CE63",  # Gold
  "Category16" = "#FABFD2",  # Light pink
  "Category17" = "#D4A6C8",  # Lilac
  "Category18" = "#D7D7D7",  # Light gray
  "Category19" = "#9D7660",  # Tan
  "Category20" = "#A0CBE8"   # Sky blue
)

debug_log("palettes loaded: primates=%d, primates_dark=%d, mammals=%d",
          length(primates_palette), length(dark_primates_palette), length(mammals_palette))
