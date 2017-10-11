# Set common formatting for all plots.
THEME_ALL <- theme(
  text=element_text(size=7),
  axis.title=element_text(size=7, face="bold"),
  axis.text=element_text(size=5.5),
  axis.line.x=element_line(color="black",size=0.5),
  axis.line.y=element_line(color="black",size=0.5),
  strip.text.x=element_text(margin = margin(3,0,3,0), size=7, face="bold"),
  strip.text.y=element_text(margin = margin(0,3,0,3), size=7, face="bold"),
  strip.background=element_blank())

# Set the standard directory for plots.
DIR <- "data/plots/"


