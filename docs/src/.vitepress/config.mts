import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import { transformerMetaWordHighlight } from '@shikijs/transformers';

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
        md.use(mathjax3),
        md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
    codeTransformers: [transformerMetaWordHighlight(),],
  },

  head: [
    //   [
    //     "script",
    //     { async: "", src: "https://www.googletagmanager.com/gtag/js?id=G-Q8GYTEVTZ2" },
    //   ],
    //   [
    //     "script",
    //     {},
    //     `window.dataLayer = window.dataLayer || [];
    //         function gtag(){dataLayer.push(arguments);}
    //         gtag('js', new Date());
    //         gtag('config', 'G-Q8GYTEVTZ2');`,
    //   ],
    ['link', { rel: 'icon', href: '/HarmonicBalance.jl/dev/favicon.ico' }],
    ['link', { rel: 'manifest', href: '/site.webmanifest' }],
  ],

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },

    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/introduction' },
      { text: 'Background', link: '/background/harmonic_balance' },
      { text: 'Tutorials', link: '/tutorials' },
      { text: 'Examples', link: '/examples' },
      {
        text: 'Manual', items: [
          { text: 'Entering equations of motion', link: '/manual/entering_eom.md' },
          { text: 'Computing effective system', link: '/manual/extracting_harmonics' },
          { text: 'Krylov-Bogoliubov', link: '/manual/Krylov-Bogoliubov_method' },
          { text: 'Time evolution', link: '/manual/time_dependent' },
          { text: 'Linear response', link: '/manual/linear_response' },
          { text: 'Plotting', link: '/manual/plotting' },
          { text: `Saving and loading`, link: '/manual/saving' },
        ]
      },
    ],

    sidebar: {
      "/introduction/": [
        { text: 'Introduction', link: '/index' },
        { text: 'Overview', link: '/introduction/overview' },
        { text: 'Tutorials', link: '/tutorials/index' },
        { text: 'Citation', link: '/introduction/citation' }
      ],
      "/background/": [
        { text: 'The method of Harmonic Balance', link: '/background/harmonic_balance.md' },
        { text: 'Stability and linear response', link: 'background/stability_response.md' },
        { text: 'Limit cycles', link: 'background/limit_cycles.md' }
      ],
      "/tutorials/": [
        { text: 'Steady states', link: '/tutorials/steady_states' },
        { text: 'Classifying solutions', link: '/tutorials/classification' },
        { text: 'Linear response', link: '/tutorials/linear_response' },
        { text: 'Transient dynamics', link: '/tutorials/time_dependent' },
        { text: 'Limit cycles', link: '/tutorials/limit_cycles' }
      ],
      "/examples/": [
        { text: 'Wave mixing', link: '/examples/wave_mixing' },
        { text: 'Parametric three wave mixing', link: '/examples/parametric_via_three_wave_mixing' },
        { text: 'Parametric oscillator', link: '/examples/parametron' }
      ],
      "/manual/": [
        { text: 'Entering equations of motion', link: '/manual/entering_eom.md' },
        { text: 'Computing effective system', link: '/manual/extracting_harmonics' },
        { text: 'Krylov-Bogoliubov', link: '/manual/Krylov-Bogoliubov_method' },
        { text: 'Time evolution', link: '/manual/time_dependent' },
        { text: 'Linear response', link: '/manual/linear_response' },
        { text: 'Plotting', link: '/manual/plotting' },
        { text: `Saving and loading`, link: '/manual/saving' }
      ],
      "/api/": []
    },

    editLink: {
      pattern: 'https://github.com/NonlinearOscillations/HarmonicBalance.jl/edit/main/docs/src/:path',
      text: 'Edit this page on GitHub'
    },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/NonlinearOscillations/HarmonicBalance.jl' },
      { icon: 'twitter', link: 'https://x.com/Zilberberg_Phys' },
    ],

    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>Released under the MIT License. Powered by the <a href="https://www.julialang.org">Julia Programming Language</a>.<br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    },
  }
})
