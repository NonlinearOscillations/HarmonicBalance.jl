import{_ as i,c as o,a5 as t,j as a,a as n,G as l,B as p,o as r}from"./chunks/framework.fPm9F4bo.js";const f=JSON.parse('{"title":"Analysis and plotting","description":"","frontmatter":{},"headers":[],"relativePath":"manual/plotting.md","filePath":"manual/plotting.md"}'),d={name:"manual/plotting.md"},c={class:"jldocstring custom-block",open:""};function h(g,s,u,k,m,y){const e=p("Badge");return r(),o("div",null,[s[3]||(s[3]=t('<h1 id="Analysis-and-plotting" tabindex="-1">Analysis and plotting <a class="header-anchor" href="#Analysis-and-plotting" aria-label="Permalink to &quot;Analysis and plotting {#Analysis-and-plotting}&quot;">​</a></h1><p>The key method for visualization is <code>transform_solutions</code>, which parses a string into a symbolic expression and evaluates it for every steady state solution.</p><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>HarmonicBalance.transform_solutions</code>. Check Documenter&#39;s build log for details.</p></div><h2 id="Plotting-solutions" tabindex="-1">Plotting solutions <a class="header-anchor" href="#Plotting-solutions" aria-label="Permalink to &quot;Plotting solutions {#Plotting-solutions}&quot;">​</a></h2><p>The function <code>plot</code> is multiple-dispatched to plot 1D and 2D datasets. In 1D, the solutions are colour-coded according to the branches obtained by <code>sort_solutions</code>.</p>',5)),a("details",c,[a("summary",null,[s[0]||(s[0]=a("a",{id:"RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}",href:"#RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}"},[a("span",{class:"jlbinding"},"RecipesBase.plot")],-1)),s[1]||(s[1]=n()),l(e,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[2]||(s[2]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result{S, P, D, F}</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    varargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cut,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p><strong>Plot a <code>Result</code> object.</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class       :   only plot solutions in this class(es) (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class   :   do not plot solutions in this class(es)</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr().</p><p>See also <code>plot!</code></p><p>The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., <code>y=2*sqrt(u1^2+v1^2)</code> plots the amplitude of the first quadratures multiplied by 2.</p><p><strong>1D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; x::String, y::String, class=&quot;default&quot;, not_class=[], kwargs...)</span></span>
<span class="line"><span>plot(res::Result, y::String; kwargs...) # take x automatically from Result</span></span></code></pre></div><p>Default behaviour is to plot stable solutions as full lines, unstable as dashed.</p><p>If a sweep in two parameters were done, i.e., <code>dim(res)==2</code>, a one dimensional cut can be plotted by using the keyword <code>cut</code> were it takes a <code>Pair{Num, Float}</code> type entry. For example, <code>plot(res, y=&quot;sqrt(u1^2+v1^2), cut=(λ =&gt; 0.2))</code> plots a cut at <code>λ = 0.2</code>.</p><hr><p><strong>2D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; z::String, branch::Int64, class=&quot;physical&quot;, not_class=[], kwargs...)</span></span></code></pre></div><p>To make the 2d plot less chaotic it is required to specify the specific <code>branch</code> to plot, labeled by a <code>Int64</code>.</p><p>The x and y axes are taken automatically from <code>res</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/ab10c6561a85e9460f9933831e9af5e83f423430/src/plotting_Plots.jl#L11" target="_blank" rel="noreferrer">source</a></p>`,17))]),s[4]||(s[4]=t('<h2 id="Plotting-phase-diagrams" tabindex="-1">Plotting phase diagrams <a class="header-anchor" href="#Plotting-phase-diagrams" aria-label="Permalink to &quot;Plotting phase diagrams {#Plotting-phase-diagrams}&quot;">​</a></h2><p>In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. <code>plot_phase_diagram</code> handles this for 1D and 2D datasets.</p><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>HarmonicBalance.plot_phase_diagram</code>. Check Documenter&#39;s build log for details.</p></div><h2 id="Plot-spaghetti-plot" tabindex="-1">Plot spaghetti plot <a class="header-anchor" href="#Plot-spaghetti-plot" aria-label="Permalink to &quot;Plot spaghetti plot {#Plot-spaghetti-plot}&quot;">​</a></h2><p>Sometimes, it is useful to plot the quadratures of the steady states (u, v) in function of a swept parameter. This is done with <code>plot_spaghetti</code>.</p><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>HarmonicBalance.plot_spaghetti</code>. Check Documenter&#39;s build log for details.</p></div>',6))])}const E=i(d,[["render",h]]);export{f as __pageData,E as default};
