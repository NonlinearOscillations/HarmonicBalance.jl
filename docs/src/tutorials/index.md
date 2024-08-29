
```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "simple_Duffing",
    src: "../assets/response_single.png",
    caption: "Steady states",
    desc: "How to get the steady states of the harmonic equations."
  },
  {
    href: "parametron",
    src: "../assets/2d_phase_diagram.png",
    caption: "Classifying solutions",
    desc: "Learn how to add different types of drives."
  },
  {
    href: "linear_response",
    src: "../assets/nonlin_F_noise.png",
    caption: "Linear response",
    desc: "Learn how to compute the linear response of a steady state."
  },
  {
    href: "time_dependent",
    src: "../assets/evo_to_steady.png",
    caption: "Stroboscopic evolution",
    desc: "Learn how to investigate stroboscopic time evolution."
  },
  {
    href: "limit_cycles",
    src: "../assets/vdp_degenerate.png",
    caption: "Limit cycles",
    desc: "Learn how to find the limit cycles of your system."
  }
];


</script>

# Tutorials

We show the capabilities of the package by providing a series of tutorials. Here we use the duffing oscillator as a typical example. Examples of other systems can be found in the example tab.

<Gallery :images="beginner" />
```