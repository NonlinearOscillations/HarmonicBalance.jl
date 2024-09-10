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
    caption: "parametric three wave mixing",
    desc: "Parametric excitation through three wave mixing."
  }
];


</script>

# Examples

<Gallery :images="examples" />
```

::: tip

If you wrote an amazing tutorial showcasing `HarmonicBalance.jl` yourself, please open an issue or PR to add it to the list!

:::