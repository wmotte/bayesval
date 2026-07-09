# Write a figure at house resolution through a device that emits plain sRGB.

R's default bitmap device on macOS (`type = "quartz"`) bakes an Apple
"Generic RGB" ICC profile into the PNG *and* converts the pixel values
into it. LaTeX ignores the embedded profile and renders the raw values,
so an Okabe-Ito orange (`#E69F00`) prints as `#DE8E08` and the blue
(`#0072B2`) as `#0B5DA2` — darker, hue-shifted, and with the palette's
contrast flattened. Pin an sRGB device so the written pixels are the
colours we asked for.

## Usage

``` r
bc_ggsave(filename, plot, width, height, dpi = 600, ...)
```

## Arguments

- filename:

  output path.

- plot:

  a ggplot object.

- width, height:

  size in inches.

- dpi:

  resolution (house style: 600).

- ...:

  passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
