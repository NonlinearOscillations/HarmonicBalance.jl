import{_ as a,c as t,a4 as i,o as e}from"./chunks/framework.1yUM_tCS.js";const g=JSON.parse('{"title":"Analysis and plotting","description":"","frontmatter":{},"headers":[],"relativePath":"manual/plotting.md","filePath":"manual/plotting.md"}'),n={name:"manual/plotting.md"};function l(o,s,p,r,d,c){return e(),t("div",null,s[0]||(s[0]=[i(`<h1 id="Analysis-and-plotting" tabindex="-1">Analysis and plotting <a class="header-anchor" href="#Analysis-and-plotting" aria-label="Permalink to &quot;Analysis and plotting {#Analysis-and-plotting}&quot;">​</a></h1><p>The key method for visualization is <code>transform_solutions</code>, which parses a string into a symbolic expression and evaluates it for every steady state solution.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="HarmonicBalance.transform_solutions" href="#HarmonicBalance.transform_solutions">#</a> <b><u>HarmonicBalance.transform_solutions</u></b> — <i>Function</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">transform_solutions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    func;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    branches,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    realify</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Vector</span></span></code></pre></div><p>Takes a <code>Result</code> object and a string <code>f</code> representing a Symbolics.jl expression. Returns an array with the values of <code>f</code> evaluated for the respective solutions. Additional substitution rules can be specified in <code>rules</code> in the format <code>(&quot;a&quot; =&gt; val)</code> or <code>(a =&gt; val)</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/33d78c7f35a4745ca425b4b8a4aceda273a5bcc5/src/transform_solutions.jl#L3" target="_blank" rel="noreferrer">source</a></p></div><br><h2 id="Plotting-solutions" tabindex="-1">Plotting solutions <a class="header-anchor" href="#Plotting-solutions" aria-label="Permalink to &quot;Plotting solutions {#Plotting-solutions}&quot;">​</a></h2><p>The function <code>plot</code> is multiple-dispatched to plot 1D and 2D datasets. In 1D, the solutions are colour-coded according to the branches obtained by <code>sort_solutions</code>.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="RecipesBase.plot-Tuple{Result, Vararg{Any}}" href="#RecipesBase.plot-Tuple{Result, Vararg{Any}}">#</a> <b><u>RecipesBase.plot</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, varargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; cut, kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p><strong>Plot a <code>Result</code> object.</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class       :   only plot solutions in this class(es) (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class   :   do not plot solutions in this class(es)</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr().</p><p>See also <code>plot!</code></p><p>The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., <code>y=2*sqrt(u1^2+v1^2)</code> plots the amplitude of the first quadratures multiplied by 2.</p><p><strong>1D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; x::String, y::String, class=&quot;default&quot;, not_class=[], kwargs...)</span></span>
<span class="line"><span>plot(res::Result, y::String; kwargs...) # take x automatically from Result</span></span></code></pre></div><p>Default behaviour is to plot stable solutions as full lines, unstable as dashed.</p><p>If a sweep in two parameters were done, i.e., <code>dim(res)==2</code>, a one dimensional cut can be plotted by using the keyword <code>cut</code> were it takes a <code>Pair{Num, Float64}</code> type entry. For example, <code>plot(res, y=&quot;sqrt(u1^2+v1^2), cut=(λ =&gt; 0.2))</code> plots a cut at <code>λ = 0.2</code>.</p><hr><p><strong>2D plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(res::Result; z::String, branch::Int64, class=&quot;physical&quot;, not_class=[], kwargs...)</span></span></code></pre></div><p>To make the 2d plot less chaotic it is required to specify the specific <code>branch</code> to plot, labeled by a <code>Int64</code>.</p><p>The x and y axes are taken automatically from <code>res</code></p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/33d78c7f35a4745ca425b4b8a4aceda273a5bcc5/src/plotting_Plots.jl#L11" target="_blank" rel="noreferrer">source</a></p></div><br><h2 id="Plotting-phase-diagrams" tabindex="-1">Plotting phase diagrams <a class="header-anchor" href="#Plotting-phase-diagrams" aria-label="Permalink to &quot;Plotting phase diagrams {#Plotting-phase-diagrams}&quot;">​</a></h2><p>In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. <code>plot_phase_diagram</code> handles this for 1D and 2D datasets.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="HarmonicBalance.plot_phase_diagram" href="#HarmonicBalance.plot_phase_diagram">#</a> <b><u>HarmonicBalance.plot_phase_diagram</u></b> — <i>Function</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Plot</span></span></code></pre></div><p>Plot the number of solutions in a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/33d78c7f35a4745ca425b4b8a4aceda273a5bcc5/src/plotting_Plots.jl#L305" target="_blank" rel="noreferrer">source</a></p></div><br><h2 id="Plot-spaghetti-plot" tabindex="-1">Plot spaghetti plot <a class="header-anchor" href="#Plot-spaghetti-plot" aria-label="Permalink to &quot;Plot spaghetti plot {#Plot-spaghetti-plot}&quot;">​</a></h2><p>Sometimes, it is useful to plot the quadratures of the steady states (u, v) in function of a swept parameter. This is done with <code>plot_spaghetti</code>.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="HarmonicBalance.plot_spaghetti" href="#HarmonicBalance.plot_spaghetti">#</a> <b><u>HarmonicBalance.plot_spaghetti</u></b> — <i>Function</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_spaghetti</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; x, y, z, kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Plot a three dimension line plot of a <code>Result</code> object as a function of the parameters. Works with 1D and 2D datasets.</p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p>Other kwargs are passed onto Plots.gr()</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/33d78c7f35a4745ca425b4b8a4aceda273a5bcc5/src/plotting_Plots.jl#L372-L384" target="_blank" rel="noreferrer">source</a></p></div><br>`,16)]))}const u=a(n,[["render",l]]);export{g as __pageData,u as default};