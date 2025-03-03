import{_ as p,C as h,c as k,o as d,ai as n,j as i,G as e,a as t,w as l}from"./chunks/framework.DHbcVXwV.js";const v=JSON.parse('{"title":"Time evolution","description":"","frontmatter":{},"headers":[],"relativePath":"manual/time_dependent.md","filePath":"manual/time_dependent.md"}'),o={name:"manual/time_dependent.md"},r={class:"jldocstring custom-block",open:""},E={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},y={class:"jldocstring custom-block",open:""};function u(m,s,b,F,f,C){const a=h("Badge");return d(),k("div",null,[s[24]||(s[24]=n('<h1 id="Time-evolution" tabindex="-1">Time evolution <a class="header-anchor" href="#Time-evolution" aria-label="Permalink to &quot;Time evolution {#Time-evolution}&quot;">​</a></h1><p>Generally, solving the ODE of oscillatory systems in time requires numerically tracking the oscillations. This is a computationally expensive process; however, using the harmonic ansatz removes the oscillatory time-dependence. Simulating instead the harmonic variables of a <code>HarmonicEquation</code> is vastly more efficient - a steady state of the system appears as a fixed point in multidimensional space rather than an oscillatory function.</p><p>The extension <code>TimeEvolution</code> is used to interface <code>HarmonicEquation</code> with the solvers contained in <code>OrdinaryDiffEq.jl</code>. Time-dependent parameter sweeps are defined using the object <code>AdiabaticSweep</code>. To use the <code>TimeEvolution</code> extension, one must first load the <code>OrdinaryDiffEq.jl</code> package.</p>',3)),i("details",r,[i("summary",null,[s[0]||(s[0]=i("a",{id:"SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent",href:"#SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent"},[i("span",{class:"jlbinding"},"SciMLBase.ODEProblem")],-1)),s[1]||(s[1]=t()),e(a,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[3]||(s[3]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ODEProblem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    eom</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    fixed_parameters;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    sweep,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    u0,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    timespan,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    perturb_initial,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Creates an ODEProblem object used by OrdinaryDiffEqTsit5.jl from the equations in <code>eom</code> to simulate time-evolution within <code>timespan</code>. <code>fixed_parameters</code> must be a dictionary mapping parameters+variables to numbers (possible to use a solution index, e.g. solutions[x][y] for branch y of solution x). If <code>u0</code> is specified, it is used as an initial condition; otherwise the values from <code>fixed_parameters</code> are used.</p>`,2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[2]||(s[2]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/ext/TimeEvolution/ODEProblem.jl#L3-L9",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),i("details",E,[i("summary",null,[s[4]||(s[4]=i("a",{id:"HarmonicBalance.AdiabaticSweep-manual-time_dependent",href:"#HarmonicBalance.AdiabaticSweep-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.AdiabaticSweep")],-1)),s[5]||(s[5]=t()),e(a,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[7]||(s[7]=n(`<p>Represents a sweep of one or more parameters of a <code>HarmonicEquation</code>. During a sweep, the selected parameters vary linearly over some timespan and are constant elsewhere.</p><p>Sweeps of different variables can be combined using <code>+</code>.</p><p><strong>Fields</strong></p><ul><li><code>functions::Dict{Num, Function}</code>: Maps each swept parameter to a function.</li></ul><p><strong>Examples</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># create a sweep of parameter a from 0 to 1 over time 0 -&gt; 100</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> @variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> a,b;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">));</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep[a](</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.5</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep[a](</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># do the same, varying two parameters simultaneously</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]], (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p>Successive sweeps can be combined,</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sweep1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.95</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2e4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sweep2 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(λ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.05</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2e4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4e4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sweep </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep1 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> sweep2</span></span></code></pre></div><p>multiple parameters can be swept simultaneously,</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sweep </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([ω </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.95</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], λ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5e-2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]], (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2e4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p>and custom sweep functions may be used.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">ωfunc</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sweep </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> AdiabaticSweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ωfunc)</span></span></code></pre></div>`,12)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[6]||(s[6]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/src/types.jl#L9-L48",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[25]||(s[25]=i("p",null,[t("In addition, one can use the "),i("code",null,"steady_state_sweep"),t(" function from "),i("code",null,"SteadyStateDiffEqExt"),t(" to perform a parameter sweep over the steady states of a system. For this one has to load "),i("code",null,"SteadyStateDiffEq.jl"),t(".")],-1)),i("details",c,[i("summary",null,[s[8]||(s[8]=i("a",{id:"HarmonicBalance.steady_state_sweep-manual-time_dependent",href:"#HarmonicBalance.steady_state_sweep-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.steady_state_sweep")],-1)),s[9]||(s[9]=t()),e(a,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[12]||(s[12]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">steady_state_sweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prob</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SteadyStateProblem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DynamicSS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; varied</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Pair</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Sweeps through a range of parameter values using a dynamic steady state solver <code>DynamicSS</code> of the <code>SteadyStateDiffEq.jl</code> package. Given a steady state problem and a parameter to vary, computes the steady state solution for each value in the sweep range. The solutions are returned as a vector where each element corresponds to the steady state found at that parameter value.</p>',2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[10]||(s[10]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/ext/SteadyStateDiffEqExt.jl#L12-L20",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1}),s[13]||(s[13]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">steady_state_sweep</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prob_np</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">NonlinearProblem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, prob_ss</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SteadyStateProblem</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">                  alg_np, alg_ss</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DynamicSS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; varied</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Pair</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Performs a parameter sweep by combining nonlinear root <code>alg_np</code> and steady state solvers <code>alg_ss</code>. For each parameter value, it first attempts a direct nonlinear root solver and checks its stability. If the solution is unstable or not found, it switches to a dynamic steady state solver. This hybrid approach is much faster then only using a steady state solver.</p>`,2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[11]||(s[11]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/ext/SteadyStateDiffEqExt.jl#L39-L47",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[26]||(s[26]=i("h2",{id:"plotting",tabindex:"-1"},[t("Plotting "),i("a",{class:"header-anchor",href:"#plotting","aria-label":'Permalink to "Plotting"'},"​")],-1)),i("details",g,[i("summary",null,[s[14]||(s[14]=i("a",{id:"RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent",href:"#RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent"},[i("span",{class:"jlbinding"},"RecipesBase.plot")],-1)),s[15]||(s[15]=t()),e(a,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[17]||(s[17]=n('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(soln</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ODESolution</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, f</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, harm_eq</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Plot a function <code>f</code> of a time-dependent solution <code>soln</code> of <code>harm_eq</code>.</p><p><strong>As a function of time</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(soln::ODESolution, f::String, harm_eq::HarmonicEquation; kwargs...)</span></span></code></pre></div><p><code>f</code> is parsed by Symbolics.jl</p><p><strong>parametric plots</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>plot(soln::ODESolution, f::Vector{String}, harm_eq::HarmonicEquation; kwargs...)</span></span></code></pre></div><p>Parametric plot of f[1] against f[2]</p><p>Also callable as plot!</p>',9)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[16]||(s[16]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/ext/PlotsExt/time_evolution.jl#L1-L19",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[27]||(s[27]=i("h2",{id:"miscellaneous",tabindex:"-1"},[t("Miscellaneous "),i("a",{class:"header-anchor",href:"#miscellaneous","aria-label":'Permalink to "Miscellaneous"'},"​")],-1)),s[28]||(s[28]=i("p",null,"Using a time-dependent simulation can verify solution stability in cases where the Jacobian is too expensive to compute.",-1)),i("details",y,[i("summary",null,[s[18]||(s[18]=i("a",{id:"HarmonicBalance.is_stable-manual-time_dependent",href:"#HarmonicBalance.is_stable-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.is_stable")],-1)),s[19]||(s[19]=t()),e(a,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[22]||(s[22]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">is_stable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    steady_solution</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">OrderedCollections.OrderedDict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    eom</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    timespan,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    tol,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    perturb_initial</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Numerically investigate the stability of a solution <code>soln</code> of <code>eom</code> within <code>timespan</code>. The initial condition is displaced by <code>perturb_initial</code>.</p><p>Return <code>true</code> the solution evolves within <code>tol</code> of the initial value (interpreted as stable).</p>`,3)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[20]||(s[20]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/ext/TimeEvolution/ODEProblem.jl#L61",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1}),s[23]||(s[23]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">is_stable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    soln</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">OrderedCollections.OrderedDict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Any</span></span></code></pre></div><p>Returns true if the solution <code>soln</code> of the Result <code>res</code> is stable. Stable solutions are real and have all Jacobian eigenvalues Re(λ) &lt;= 0. <code>im_tol</code> : an absolute threshold to distinguish real/complex numbers. <code>rel_tol</code>: Re(λ) considered &lt;=0 if real.(λ) &lt; rel_tol*abs(λmax)</p>`,2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[21]||(s[21]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/4be0f9d09f27c00c7d434acffd8dd8d71fa6448d/src/classification.jl#L77",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})])])}const A=p(o,[["render",u]]);export{v as __pageData,A as default};
