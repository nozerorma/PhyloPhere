# Palettes

# Comment, uncomment and add families as needed depending on the scope of the analysis.

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

## Desaturated primate color palette
primates_palette <- c(
  "Aotidae"           = "#6BB9B3",
  "Callitrichidae"    = "#6F94B5",
  "Cebidae"           = "#D95F5C",
  "Atelidae"          = "#B0C7D2",
  "Cercopithecidae"   = "#70A774",
  "Hominidae"         = "#F4B0B0",
  "Hylobatidae"       = "#7C6FAA",
  "Lemuridae"         = "#F89F4E",
  "Indriidae"         = "#C6DFBC",
  "Cheirogaleidae"    = "#F7C998",
  "Lorisidae"         = "#A78CB3",
  "Galagidae"         = "#E0C769",
  "Pitheciidae"       = "#B57E4B",
  "Daubentoniidae"    = "#A65D8A",
  "Tarsiidae"         = "#C3A6D9"
)

## Darkened primate color palette
dark_primates_palette <- c(
  "Aotidae"         = "#1B5C58",
  "Callitrichidae"  = "#2F495E",
  "Cebidae"         = "#7A221F",
  "Atelidae"        = "#41626E",
  "Cercopithecidae" = "#305333",
  "Hominidae"       = "#AD2A2A",
  "Hylobatidae"     = "#3E365B",
  "Lemuridae"       = "#7F4900",
  "Indriidae"       = "#497033",
  "Cheirogaleidae"  = "#865B00",
  "Lorisidae"       = "#56415F",
  "Galagidae"       = "#6F5E03",
  "Pitheciidae"     = "#5C3D20",
  "Daubentoniidae"  = "#572C47",
  "Tarsiidae"       = "#6D4485"
)

# Desaturated mammalian color palette
mammals_palette <- c(
  "Carnivora"        = "#4B6F87",
  "Perissodactyla"   = "#A65657",
  "Hyracoidea"       = "#B8C8D0",
  "Afrosoricida"     = "#517D4E",
  "Artiodactyla"     = "#DDB6B6",
  "Rodentia"         = "#D2B999",
  "Chiroptera"       = "#B27F4C",
  "Diprotodontia"    = "#B3C5A3",
  "Scandentia"       = "#C6BCCB",
  "Proboscidea"      = "#9CB07B",
  "Sirenia"          = "#986742",
  "Tubulidentata"    = "#86849E",
  "Macroscelidea"    = "#AE6188",
  "Eulipotyphla"     = "#637D46",
  "Lagomorpha"       = "#A18A46",
  "Pholidota"        = "#7C6946",
  "Dermoptera"       = "#666666"
)

## Compatibility aliases
primate_family_colors <- primates_palette
dark_family_palette <- dark_primates_palette
mammal_order_colors <- mammals_palette

## Include below as needed
# Fallback palette
fallback_palette <- c(
  "Category1" = "#1f77b4",
  "Category2" = "#ff7f0e",
  "Category3" = "#2ca02c",
  "Category4" = "#d62728",
  "Category5" = "#9467bd",
  "Category6" = "#8c564b",
  "Category7" = "#e377c2",
  "Category8" = "#7f7f7f",
  "Category9" = "#bcbd22",
  "Category10"= "#17becf",
  "Category11"= "#aec7e8",
  "Category12"= "#ffbb78",
  "Category13"= "#98df8a",
  "Category14"= "#ff9896",
  "Category15"= "#c5b0d5",
  "Category16"= "#c49c94",
  "Category17"= "#f7b6d2",
  "Category18"= "#c7c7c7",
  "Category19"= "#dbdb8d",
  "Category20"= "#9edae5"
)

debug_log("palettes loaded: primates=%d, primates_dark=%d, mammals=%d",
          length(primates_palette), length(dark_primates_palette), length(mammals_palette))
