import{_ as t,C as l,c,o as d,ai as n,j as e,a as i,G as o}from"./chunks/framework.CFSb-KsM.js";const k=JSON.parse('{"title":"Saving and loading","description":"","frontmatter":{},"headers":[],"relativePath":"manual/saving.md","filePath":"manual/saving.md"}'),r={name:"manual/saving.md"},p={class:"jldocstring custom-block",open:""},m={class:"jldocstring custom-block",open:""},u={class:"jldocstring custom-block",open:""};function g(h,a,b,v,_,f){const s=l("Badge");return d(),c("div",null,[a[9]||(a[9]=n('<h1 id="Saving-and-loading" tabindex="-1">Saving and loading <a class="header-anchor" href="#Saving-and-loading" aria-label="Permalink to &quot;Saving and loading {#Saving-and-loading}&quot;">​</a></h1><p>All of the types native to <code>HarmonicBalance.jl</code> can be saved into a <code>.jld2</code> file using <code>save</code> and loaded using <code>load</code>. Most of the saving/loading is performed using the package <code>JLD2.jl</code>, with the addition of reinstating the symbolic variables in the <code>HarmonicBalance</code> namespace (needed to parse expressions used in the plotting functions) and recompiling stored functions (needed to evaluate Jacobians). As a consequence, composite objects such as <a href="/HarmonicBalance.jl/previews/PR366/manual/solving_harmonics#HarmonicBalance.Result-manual-solving_harmonics"><code>HarmonicBalance.Result</code></a> can be saved and loaded with no loss of information.</p><p>The function <code>export_csv</code> saves a .csv file which can be plot elsewhere.</p>',3)),e("details",p,[e("summary",null,[a[0]||(a[0]=e("a",{id:"HarmonicBalance.save-manual-saving",href:"#HarmonicBalance.save-manual-saving"},[e("span",{class:"jlbinding"},"HarmonicBalance.save")],-1)),a[1]||(a[1]=i()),o(s,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),a[2]||(a[2]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">save</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(filename, object)</span></span></code></pre></div><p>Saves <code>object</code> into <code>.jld2</code> file <code>filename</code> (the suffix is added automatically if not entered). The resulting file contains a dictionary with a single entry.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/957701923ba2acf71be34ea95fdd14e79b605f2a/src/saving.jl#L1-L7" target="_blank" rel="noreferrer">source</a></p>',3))]),e("details",m,[e("summary",null,[a[3]||(a[3]=e("a",{id:"HarmonicBalance.load-manual-saving",href:"#HarmonicBalance.load-manual-saving"},[e("span",{class:"jlbinding"},"HarmonicBalance.load")],-1)),a[4]||(a[4]=i()),o(s,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),a[5]||(a[5]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">load</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(filename)</span></span></code></pre></div><p>Loads an object from <code>filename</code>. For objects containing symbolic expressions such as <code>HarmonicEquation</code>, the symbolic variables are reinstated in the <code>HarmonicBalance</code> namespace.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/957701923ba2acf71be34ea95fdd14e79b605f2a/src/saving.jl#L22-L28" target="_blank" rel="noreferrer">source</a></p>',3))]),e("details",u,[e("summary",null,[a[6]||(a[6]=e("a",{id:"HarmonicBalance.export_csv-manual-saving",href:"#HarmonicBalance.export_csv-manual-saving"},[e("span",{class:"jlbinding"},"HarmonicBalance.export_csv")],-1)),a[7]||(a[7]=i()),o(s,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),a[8]||(a[8]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">export_csv</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(filename, res, branch)</span></span></code></pre></div><p>Saves into <code>filename</code> a specified solution <code>branch</code> of the Result <code>res</code>.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/957701923ba2acf71be34ea95fdd14e79b605f2a/src/saving.jl#L77-L81" target="_blank" rel="noreferrer">source</a></p>',3))])])}const T=t(r,[["render",g]]);export{k as __pageData,T as default};
