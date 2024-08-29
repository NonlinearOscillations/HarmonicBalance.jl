
```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "tutorials/simple_Duffing",
    src: "../assets/response_single.png",
    caption: "Finding the steady states",
    desc: "How to get started with Julia and HarmonicBalance for those who have never used Julia before."
  },
  {
    href: "tutorials/linear_response",
    src: "../assets/nonlin_F_noise.png",
    caption: "Compute the linear response",
    desc: "Learn how to compute the linear response of a steady state."
  },
  {
    href: "tutorials/time_dependent",
    src: "../assets/evo_to_steady.png",
    caption: "Evolve the stroboscopic system in time",
    desc: "Learn how to investigate custom layers and train an RNN on time-series data."
  },
  {
    href: "tutorials/limit_cycles",
    src: "../assets/vdp_degenerate.png",
    caption: "Find limit cycles",
    desc: "Learn how to find the limit cycles of your system."
  },
  {
    href: "tutorials/parametron",
    src: "../assets/2d_phase_diagram.png",
    caption: "Adding parametric driving",
    desc: "Learn how to add different types of drives."
  }
];


</script>

# Tutorials

We show the capabilities of the package by providing a series of tutorials. Here we use the duffing oscillator as a typical example. Examples of other systems can be found in the example tab.

<Gallery :images="beginner" />
```