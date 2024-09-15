```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const examples = [
  {
    href: "wave_mixing",
    src: "",
    caption: "Wave mixing",
    desc: "Understand three and four wave mixing."
  },
  {
    href: "parametric_via_three_wave_mixing",
    src: "",
    caption: "Parametric three wave mixing",
    desc: "Parametric excitation through three wave mixing."
  },
  {
    href: "parametron",
    src: "",
    caption: "Parametric oscillator",
    desc: "Introduction to the parametric oscillator."
  }
];


</script>
```
# [Examples](@id examples)

```@raw html
<Gallery :images="examples" />
```

::: tip

If you wrote an amazing tutorial showcasing `HarmonicBalance.jl` yourself, please open an issue or PR to add it to the list!

:::