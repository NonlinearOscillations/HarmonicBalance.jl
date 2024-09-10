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
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
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
      "/introduction/": {
        text: 'Getting Started', collapsed: false, items: [
          { text: 'Introduction', link: '/index' },
          { text: 'Overview', link: '/introduction/overview' },
          { text: 'Resources', link: '/introduction/resources' },
          { text: 'Citation', link: '/introduction/citation' }]
      },
      "/background/": {
        text: 'Background', collapsed: false, items: [
          { text: 'The method of Harmonic Balance', link: '/background/harmonic_balance.md' },
          { text: 'Stability and linear response', link: 'background/stability_response.md' },
          { text: 'Limit cycles', link: 'background/limit_cycles.md' },
        ]
      },
      "/tutorials/": {
        text: 'Tutorials', collapsed: false, items: [
          { text: 'Steady states', link: '/tutorials/duffing' },
          { text: 'Transient dynamics', link: '/tutorials/time_dependent' },
          { text: 'Classifying solutions', link: '/tutorials/parametron' },
          { text: 'Linear response', link: '/tutorials/linear_response' },
          { text: 'Limit cycle', link: '/tutorials/limit_cycles' },
        ]
      },
      "/examples/": {
        text: 'Examples', collapsed: false, items: [
          { text: 'Wave mixing', link: '/examples/wave_mixing' },
        ]
      },
      "/manual/": {
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
      "/api/": {
        text: 'API Reference', collapsed: false, items: []
      }
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
