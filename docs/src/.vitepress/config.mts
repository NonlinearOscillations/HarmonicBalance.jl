import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import { transformerMetaWordHighlight } from '@shikijs/transformers';

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // TODO: replace this in makedocs!
  title: 'HarmonicBalance.jl',
  description: 'A Julia package for solving nonlinear differential equations using the harmonic balance method.',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: '/HarmonicBalance.jl/dev/favicon.ico' }]],

  vite: {
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    }
  },

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
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    logo: { src: '/logo.png', width: 24, height: 24 },

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
      {
        text: 'Tutorials', items: [
          { text: 'Steady states', link: '/examples/simple_Duffing.md' },
          { text: 'Transient dynamics', link: '/examples/time_dependent.md' },
          { text: 'Classifying solutions', link: '/examples/parametron.md' },
          { text: 'Linear response', link: '/examples/linear_response.md' },
          { text: 'Limit cycle', link: '/examples/limit_cycles.md' },
        ]
      },
      { text: 'Examples', link: '/examples/overview' },
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
          { text: 'Steady states', link: '/examples/simple_Duffing' },
          { text: 'Transient dynamics', link: '/examples/time_dependent' },
          { text: 'Classifying solutions', link: '/examples/parametron' },
          { text: 'Linear response', link: '/examples/linear_response' },
          { text: 'Limit cycle', link: '/examples/limit_cycles' },
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
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})