import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import { transformerMetaWordHighlight } from '@shikijs/transformers';

const baseTemp = {
  base: '/HarmonicBalance.jl/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', link: '/introduction/index' },
{ text: 'Theoretical Background', collapsed: false, items: [
{ text: 'Floquet expansions', link: '/background/harmonic_balance' },
{ text: 'Stability and Linear Response', link: '/background/stability_response' },
{ text: 'Limit Cycles', link: '/background/limit_cycles' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Steady States', link: '/tutorials/steady_states' },
{ text: 'Classifying Solutions', link: '/tutorials/classification' },
{ text: 'Linear Response', link: '/tutorials/linear_response' },
{ text: 'Transient Dynamics', link: '/tutorials/time_dependent' },
{ text: 'Limit Cycles', link: '/tutorials/limit_cycles' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Wave Mixing', link: '/examples/wave_mixing' },
{ text: 'Parametric Three Wave Mixing', link: '/examples/parametric_via_three_wave_mixing' },
{ text: 'Parametric Oscillator', link: '/examples/parametron' },
{ text: 'State Dependent Perturbation', link: '/examples/state_dependent_perturbation' }]
 }]
 },
{ text: 'Resources', collapsed: false, items: [
{ text: 'API', link: '/manual/API' },
{ text: 'Bibliography', link: '/introduction/citation' },
{ text: 'Manual', collapsed: false, items: [
{ text: 'Entering Equations of Motion', link: '/manual/entering_eom' },
{ text: 'Computing Effective System', link: '/manual/extracting_harmonics' },
{ text: 'Computing Steady States', link: '/manual/methods' },
{ text: 'Krylov-Bogoliubov', link: '/manual/Krylov-Bogoliubov_method' },
{ text: 'Time Evolution', link: '/manual/time_dependent' },
{ text: 'Linear Response', link: '/manual/linear_response' },
{ text: 'Plotting', link: '/manual/plotting' },
{ text: 'Saving and Loading', link: '/manual/saving' }]
 }]
 }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: baseTemp.base,
  title: 'HarmonicBalance.jl',
  description: 'Documentation for HarmonicBalance.jl',
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [
    [
      "script",
      { async: "", src: "https://www.googletagmanager.com/gtag/js?id=G-RE962QZ6DQ" },
    ],
    [
      "script",
      {},
      `window.dataLayer = window.dataLayer || [];
              function gtag(){dataLayer.push(arguments);}
              gtag('js', new Date());
              gtag('config', 'G-RE962QZ6DQ');`,
    ],
    ['link', { rel: 'icon', href: '/HarmonicBalance.jl/dev/favicon.ico' }],
    ['link', { rel: 'icon', href: '/favicon.ico' }],
    ['link', { rel: 'manifest', href: '/site.webmanifest' }],

    ['script', { src: `/HarmonicBalance.jl/versions.js` }],
    ['script', { src: `${baseTemp.base}siteinfo.js` }]
  ],
  ignoreDeadLinks: true,

  markdown: {
    math: true,

    // options for @mdit-vue/plugin-toc
    // https://github.com/mdit-vue/mdit-vue/tree/main/packages/plugin-toc#options
    toc: { level: [2, 3, 4] }, // for API page, triggered by: [[toc]]

    config(md) {
      md.use(tabsMarkdownPlugin),
        md.use(mathjax3),
        md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    }
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', link: '/introduction/index' },
{ text: 'Theoretical Background', collapsed: false, items: [
{ text: 'Floquet expansions', link: '/background/harmonic_balance' },
{ text: 'Stability and Linear Response', link: '/background/stability_response' },
{ text: 'Limit Cycles', link: '/background/limit_cycles' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Steady States', link: '/tutorials/steady_states' },
{ text: 'Classifying Solutions', link: '/tutorials/classification' },
{ text: 'Linear Response', link: '/tutorials/linear_response' },
{ text: 'Transient Dynamics', link: '/tutorials/time_dependent' },
{ text: 'Limit Cycles', link: '/tutorials/limit_cycles' },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Wave Mixing', link: '/examples/wave_mixing' },
{ text: 'Parametric Three Wave Mixing', link: '/examples/parametric_via_three_wave_mixing' },
{ text: 'Parametric Oscillator', link: '/examples/parametron' },
{ text: 'State Dependent Perturbation', link: '/examples/state_dependent_perturbation' }]
 }]
 },
{ text: 'Resources', collapsed: false, items: [
{ text: 'API', link: '/manual/API' },
{ text: 'Bibliography', link: '/introduction/citation' },
{ text: 'Manual', collapsed: false, items: [
{ text: 'Entering Equations of Motion', link: '/manual/entering_eom' },
{ text: 'Computing Effective System', link: '/manual/extracting_harmonics' },
{ text: 'Computing Steady States', link: '/manual/methods' },
{ text: 'Krylov-Bogoliubov', link: '/manual/Krylov-Bogoliubov_method' },
{ text: 'Time Evolution', link: '/manual/time_dependent' },
{ text: 'Linear Response', link: '/manual/linear_response' },
{ text: 'Plotting', link: '/manual/plotting' },
{ text: 'Saving and Loading', link: '/manual/saving' }]
 }]
 }
]
,
    editLink: { pattern: "https://github.com/NonlinearOscillations/HarmonicBalance.jl/edit/master/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/NonlinearOscillations/HarmonicBalance.jl' },
      { icon: 'twitter', link: 'https://x.com/Zilberberg_Phys' },
    ],

    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    },
  }
})