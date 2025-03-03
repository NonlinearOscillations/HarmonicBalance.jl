import{_ as l,c as o,j as a,a as t,G as e,a4 as n,B as p,o as r}from"./chunks/framework.Ask75_is.js";const j=JSON.parse('{"title":"Analysis and plotting","description":"","frontmatter":{},"headers":[],"relativePath":"manual/plotting.md","filePath":"manual/plotting.md"}'),d={name:"manual/plotting.md"},h={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""};function u(b,s,y,m,E,f){const i=p("Badge");return r(),o("div",null,[s[12]||(s[12]=a("h1",{id:"Analysis-and-plotting",tabindex:"-1"},[t("Analysis and plotting "),a("a",{class:"header-anchor",href:"#Analysis-and-plotting","aria-label":'Permalink to "Analysis and plotting {#Analysis-and-plotting}"'},"​")],-1)),s[13]||(s[13]=a("p",null,[t("The key method for visualization is "),a("code",null,"transform_solutions"),t(", which parses a string into a symbolic expression and evaluates it for every steady state solution.")],-1)),a("details",h,[a("summary",null,[s[0]||(s[0]=a("a",{id:"HarmonicBalance.transform_solutions",href:"#HarmonicBalance.transform_solutions"},[a("span",{class:"jlbinding"},"HarmonicBalance.transform_solutions")],-1)),s[1]||(s[1]=t()),e(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[2]||(s[2]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">transform_solutions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    func;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    branches,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    realify</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Vector</span></span></code></pre></div><p>Takes a <code>Result</code> object and a string <code>f</code> representing a Symbolics.jl expression. Returns an array with the values of <code>f</code> evaluated for the respective solutions. Additional substitution rules can be specified in <code>rules</code> in the format <code>(&quot;a&quot; =&gt; val)</code> or <code>(a =&gt; val)</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/896553e4bbbd87df4a94995c58e14bd9d25d4976/src/transform_solutions.jl#L38" target="_blank" rel="noreferrer">source</a></p>`,3))]),s[14]||(s[14]=a("h2",{id:"Plotting-solutions",tabindex:"-1"},[t("Plotting solutions "),a("a",{class:"header-anchor",href:"#Plotting-solutions","aria-label":'Permalink to "Plotting solutions {#Plotting-solutions}"'},"​")],-1)),s[15]||(s[15]=a("p",null,[t("The function "),a("code",null,"plot"),t(" is multiple-dispatched to plot 1D and 2D datasets. In 1D, the solutions are colour-coded according to the branches obtained by "),a("code",null,"sort_solutions"),t(".")],-1)),a("details",c,[a("summary",null,[s[3]||(s[3]=a("a",{id:"RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}",href:"#RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}"},[a("span",{class:"jlbinding"},"RecipesBase.plot")],-1)),s[4]||(s[4]=t()),e(i,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[5]||(s[5]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    varargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cut,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p><strong>Plot a <code>Result</code> object.</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class       :   only plot solutions in this class(es) (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class   :   do not plot solutions in this class(es)</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr().</p><p>See also <code>plot!</code></p><p>The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., <code>y=2*sqrt(u1^2+v1^2)</code> plots the amplitude of the first quadratures multiplied by 2.</p><p><strong>1D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; x::String, y::String, class=&quot;default&quot;, not_class=[], kwargs...)</span></span>
<span class="line"><span>plot(res::Result, y::String; kwargs...) # take x automatically from Result</span></span></code></pre></div><p>Default behaviour is to plot stable solutions as full lines, unstable as dashed.</p><p>If a sweep in two parameters were done, i.e., <code>dim(res)==2</code>, a one dimensional cut can be plotted by using the keyword <code>cut</code> were it takes a <code>Pair{Num, Float64}</code> type entry. For example, <code>plot(res, y=&quot;sqrt(u1^2+v1^2), cut=(λ =&gt; 0.2))</code> plots a cut at <code>λ = 0.2</code>.</p><hr><p><strong>2D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; z::String, branch::Int64, class=&quot;physical&quot;, not_class=[], kwargs...)</span></span></code></pre></div><p>To make the 2d plot less chaotic it is required to specify the specific <code>branch</code> to plot, labeled by a <code>Int64</code>.</p><p>The x and y axes are taken automatically from <code>res</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/896553e4bbbd87df4a94995c58e14bd9d25d4976/src/plotting_Plots.jl#L11" target="_blank" rel="noreferrer">source</a></p>`,17))]),s[16]||(s[16]=a("h2",{id:"Plotting-phase-diagrams",tabindex:"-1"},[t("Plotting phase diagrams "),a("a",{class:"header-anchor",href:"#Plotting-phase-diagrams","aria-label":'Permalink to "Plotting phase diagrams {#Plotting-phase-diagrams}"'},"​")],-1)),s[17]||(s[17]=a("p",null,[t("In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. "),a("code",null,"plot_phase_diagram"),t(" handles this for 1D and 2D datasets.")],-1)),a("details",g,[a("summary",null,[s[6]||(s[6]=a("a",{id:"HarmonicBalance.plot_phase_diagram",href:"#HarmonicBalance.plot_phase_diagram"},[a("span",{class:"jlbinding"},"HarmonicBalance.plot_phase_diagram")],-1)),s[7]||(s[7]=t()),e(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[8]||(s[8]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p>Plot the number of solutions in a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/896553e4bbbd87df4a94995c58e14bd9d25d4976/src/plotting_Plots.jl#L274" target="_blank" rel="noreferrer">source</a></p>`,6))]),s[18]||(s[18]=a("h2",{id:"Plot-spaghetti-plot",tabindex:"-1"},[t("Plot spaghetti plot "),a("a",{class:"header-anchor",href:"#Plot-spaghetti-plot","aria-label":'Permalink to "Plot spaghetti plot {#Plot-spaghetti-plot}"'},"​")],-1)),s[19]||(s[19]=a("p",null,[t("Sometimes, it is useful to plot the quadratures of the steady states (u, v) in function of a swept parameter. This is done with "),a("code",null,"plot_spaghetti"),t(".")],-1)),a("details",k,[a("summary",null,[s[9]||(s[9]=a("a",{id:"HarmonicBalance.plot_spaghetti",href:"#HarmonicBalance.plot_spaghetti"},[a("span",{class:"jlbinding"},"HarmonicBalance.plot_spaghetti")],-1)),s[10]||(s[10]=t()),e(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[11]||(s[11]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_spaghetti</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; x, y, z, kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Plot a three dimension line plot of a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/896553e4bbbd87df4a94995c58e14bd9d25d4976/src/plotting_Plots.jl#L341-L353" target="_blank" rel="noreferrer">source</a></p>`,6))])])}const C=l(d,[["render",u]]);export{j as __pageData,C as default};
