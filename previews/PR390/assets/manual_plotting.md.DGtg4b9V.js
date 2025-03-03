import{_ as l,C as p,c as o,o as h,ai as t,j as a,a as i,G as e}from"./chunks/framework.B4T3iM86.js";const v=JSON.parse('{"title":"Plotting solutions","description":"","frontmatter":{},"headers":[],"relativePath":"manual/plotting.md","filePath":"manual/plotting.md"}'),r={name:"manual/plotting.md"},d={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""};function u(E,s,y,m,b,f){const n=p("Badge");return h(),o("div",null,[s[12]||(s[12]=t('<h1 id="plotting" tabindex="-1">Plotting solutions <a class="header-anchor" href="#plotting" aria-label="Permalink to &quot;Plotting solutions {#plotting}&quot;">​</a></h1><p>HarmonicBalance.jl comes with a plotting module <code>PlotsExt</code> that allows you to visualize the steady states in the <a href="/HarmonicBalance.jl/previews/PR390/manual/solving_harmonics#HarmonicBalance.Result-manual-solving_harmonics"><code>HarmonicBalance.Result</code></a>. The module is conditionally loaded based on the <code>Plots.jl</code> package being loaded.</p><p>The function <code>plot</code> is multiple-dispatched to plot 1D and 2D datasets. In 1D, the solutions are colour-coded according to the branches obtained by <code>sort_solutions</code>.</p>',3)),a("details",d,[a("summary",null,[s[0]||(s[0]=a("a",{id:"RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}-manual-plotting",href:"#RecipesBase.plot-Tuple{HarmonicBalance.Result, Vararg{Any}}-manual-plotting"},[a("span",{class:"jlbinding"},"RecipesBase.plot")],-1)),s[1]||(s[1]=i()),e(n,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[2]||(s[2]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result{D, S, P, F}</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    varargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    cut,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p><strong>Plot a <code>Result</code> object.</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class       :   only plot solutions in this class(es) (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class   :   do not plot solutions in this class(es)</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr().</p><p>See also <code>plot!</code></p><p>The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., <code>y=2*sqrt(u1^2+v1^2)</code> plots the amplitude of the first quadratures multiplied by 2.</p><p><strong>1D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; x::String, y::String, class=&quot;default&quot;, not_class=[], kwargs...)</span></span>
<span class="line"><span>plot(res::Result, y::String; kwargs...) # take x automatically from Result</span></span></code></pre></div><p>Default behaviour is to plot stable solutions as full lines, unstable as dashed.</p><p>If a sweep in two parameters were done, i.e., <code>dimension(res)==2</code>, a one dimensional cut can be plotted by using the keyword <code>cut</code> were it takes a <code>Pair{Num, Float}</code> type entry. For example, <code>plot(res, y=&quot;sqrt(u1^2+v1^2), cut=(λ =&gt; 0.2))</code> plots a cut at <code>λ = 0.2</code>.</p><p><strong>2D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; z::String, branch::Int64, class=&quot;physical&quot;, not_class=[], kwargs...)</span></span></code></pre></div><p>To make the 2d plot less chaotic it is required to specify the specific <code>branch</code> to plot, labeled by a <code>Int64</code>.</p><p>The x and y axes are taken automatically from <code>res</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/aa0fbc9b783f3318f8f484165d90fb3dec0c720f/ext/PlotsExt/steady_states.jl#L2" target="_blank" rel="noreferrer">source</a></p>`,16))]),a("details",c,[a("summary",null,[s[3]||(s[3]=a("a",{id:"RecipesBase.plot!-manual-plotting",href:"#RecipesBase.plot!-manual-plotting"},[a("span",{class:"jlbinding"},"RecipesBase.plot!")],-1)),s[4]||(s[4]=i()),e(n,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[5]||(s[5]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    varargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p>Similar to <code>plot</code> but adds a plot onto an existing plot.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/aa0fbc9b783f3318f8f484165d90fb3dec0c720f/ext/PlotsExt/steady_states.jl#L53" target="_blank" rel="noreferrer">source</a></p>`,3))]),s[13]||(s[13]=a("h2",{id:"Plotting-phase-diagrams",tabindex:"-1"},[i("Plotting phase diagrams "),a("a",{class:"header-anchor",href:"#Plotting-phase-diagrams","aria-label":'Permalink to "Plotting phase diagrams {#Plotting-phase-diagrams}"'},"​")],-1)),s[14]||(s[14]=a("p",null,[i("In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. "),a("code",null,"plot_phase_diagram"),i(" handles this for 1D and 2D datasets.")],-1)),a("details",k,[a("summary",null,[s[6]||(s[6]=a("a",{id:"HarmonicBalance.plot_phase_diagram-manual-plotting",href:"#HarmonicBalance.plot_phase_diagram-manual-plotting"},[a("span",{class:"jlbinding"},"HarmonicBalance.plot_phase_diagram")],-1)),s[7]||(s[7]=i()),e(n,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[8]||(s[8]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result{D, SolType}</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> SolType</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Number</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p>Plot the number of solutions in a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/aa0fbc9b783f3318f8f484165d90fb3dec0c720f/ext/PlotsExt/steady_states.jl#L231" target="_blank" rel="noreferrer">source</a></p>`,6))]),s[15]||(s[15]=a("h2",{id:"Plot-spaghetti-plot",tabindex:"-1"},[i("Plot spaghetti plot "),a("a",{class:"header-anchor",href:"#Plot-spaghetti-plot","aria-label":'Permalink to "Plot spaghetti plot {#Plot-spaghetti-plot}"'},"​")],-1)),s[16]||(s[16]=a("p",null,[i("Sometimes, it is useful to plot the quadratures of the steady states (u, v) in function of a swept parameter. This is done with "),a("code",null,"plot_spaghetti"),i(".")],-1)),a("details",g,[a("summary",null,[s[9]||(s[9]=a("a",{id:"HarmonicBalance.plot_spaghetti-manual-plotting",href:"#HarmonicBalance.plot_spaghetti-manual-plotting"},[a("span",{class:"jlbinding"},"HarmonicBalance.plot_spaghetti")],-1)),s[10]||(s[10]=i()),e(n,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[11]||(s[11]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_spaghetti</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result{D, SolType}</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> SolType</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Number</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    y,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    z,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    class,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    not_class,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    add,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Plot a three dimension line plot of a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/aa0fbc9b783f3318f8f484165d90fb3dec0c720f/ext/PlotsExt/steady_states.jl#L298" target="_blank" rel="noreferrer">source</a></p>`,6))])])}const F=l(r,[["render",u]]);export{v as __pageData,F as default};
