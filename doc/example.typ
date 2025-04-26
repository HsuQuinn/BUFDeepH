#import "lib.typ": *
#show: template.with(
  title: [Siesta Band Unfolding Document],
  short_title: "manual",
  description: [
    A manual for the BUFDeepH code.
  ],
  authors: (
    (
      name: "Quinn Hsu",
      link: "https://github.com/HsuQuinn",
    ),
  ),
  paper_size: "a4",
  // landscape: true,
  cols: 1,
  
  text_font: "Proxima Nova",
  code_font: "Cascadia Mono",
  // accent: "#3666FA", // blue
  accent: "#A68219", // blue
  h1-prefix: "Section",
  colortab: true,
)

#include "content/doc.typ"
