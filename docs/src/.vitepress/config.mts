import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import { transformerMetaWordHighlight } from '@shikijs/transformers';

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
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
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
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
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
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
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl' },
      { icon: 'twitter', link: 'https://x.com/Zilberberg_Phys' },
    ],

    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    },
  }
})