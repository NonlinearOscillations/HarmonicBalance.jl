import{_ as e,c as h,j as s,a,G as t,a4 as l,B as k,o as p}from"./chunks/framework.bmw6wTPo.js";const b=JSON.parse('{"title":"Entering equations of motion","description":"","frontmatter":{},"headers":[],"relativePath":"manual/entering_eom.md","filePath":"manual/entering_eom.md"}'),r={name:"manual/entering_eom.md"},d={class:"jldocstring custom-block",open:""},E={class:"jldocstring custom-block",open:""},o={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""};function y(c,i,F,u,f,m){const n=k("Badge");return p(),h("div",null,[i[12]||(i[12]=s("h1",{id:"Entering-equations-of-motion",tabindex:"-1"},[a("Entering equations of motion "),s("a",{class:"header-anchor",href:"#Entering-equations-of-motion","aria-label":'Permalink to "Entering equations of motion {#Entering-equations-of-motion}"'},"​")],-1)),i[13]||(i[13]=s("p",null,[a("The struct "),s("code",null,"DifferentialEquation"),a(" is the primary input method; it holds an ODE or a coupled system of ODEs composed of terms with harmonic time-dependence The dependent variables are specified during input, any other symbols are identified as parameters. Information on which variable is to be expanded in which harmonic is specified using "),s("code",null,"add_harmonic!"),a(".")],-1)),i[14]||(i[14]=s("p",null,[s("code",null,"DifferentialEquation.equations"),a(" stores a dictionary assigning variables to equations. This information is necessary because the harmonics belonging to a variable are later used to Fourier-transform its corresponding ODE.")],-1)),s("details",d,[s("summary",null,[i[0]||(i[0]=s("a",{id:"HarmonicBalance.DifferentialEquation",href:"#HarmonicBalance.DifferentialEquation"},[s("span",{class:"jlbinding"},"HarmonicBalance.DifferentialEquation")],-1)),i[1]||(i[1]=a()),t(n,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[2]||(i[2]=l(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">mutable struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DifferentialEquation</span></span></code></pre></div><p>Holds differential equation(s) of motion and a set of harmonics to expand each variable. This is the primary input for <code>HarmonicBalance.jl</code> ; after inputting the equations, the harmonics ansatz needs to be specified using <code>add_harmonic!</code>.</p><p><strong>Fields</strong></p><ul><li><p><code>equations::OrderedCollections.OrderedDict{Num, Equation}</code>: Assigns to each variable an equation of motion.</p></li><li><p><code>harmonics::OrderedCollections.OrderedDict{Num, OrderedCollections.OrderedSet{Num}}</code>: Assigns to each variable a set of harmonics.</p></li></ul><p><strong>Example</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> @variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> t, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">y</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t), ω0, ω, F, k;</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># equivalent ways to enter the simple harmonic oscillator</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x,t,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">t), x);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x,t,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">t), x);</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># two coupled oscillators, one of them driven</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x,t,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> k</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">y, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(y,t,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> y </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> k</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">x] </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [F </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">t), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], [x,y]);</span></span></code></pre></div><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/afa425ceb096439a4ff601e5145532d557edf0e3/src/types.jl#L7" target="_blank" rel="noreferrer">source</a></p>`,7))]),s("details",E,[s("summary",null,[i[3]||(i[3]=s("a",{id:"HarmonicBalance.add_harmonic!",href:"#HarmonicBalance.add_harmonic!"},[s("span",{class:"jlbinding"},"HarmonicBalance.add_harmonic!")],-1)),i[4]||(i[4]=a()),t(n,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[5]||(i[5]=l(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add_harmonic!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(diff_eom</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, var</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Num</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, ω)</span></span></code></pre></div><p>Add the harmonic <code>ω</code> to the harmonic ansatz used to expand the variable <code>var</code> in <code>diff_eom</code>.</p><p><strong>Example</strong></p><p><strong>define the simple harmonic oscillator and specify that x(t) oscillates with frequency ω</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> @variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> t, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">y</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t), ω0, ω, F, k;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> diff_eq </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x,t,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">t), x);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> add_harmonic!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(diff_eq, x, ω) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># expand x using ω</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">System of </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> differential equations</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Variables</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">       x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Harmonic ansatz</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω;</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ω0</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Differential</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t)(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Differential</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t)(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t))) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ω)</span></span></code></pre></div><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/afa425ceb096439a4ff601e5145532d557edf0e3/src/DifferentialEquation.jl#L1" target="_blank" rel="noreferrer">source</a></p>`,6))]),s("details",o,[s("summary",null,[i[6]||(i[6]=s("a",{id:"Symbolics.get_variables-Tuple{DifferentialEquation}",href:"#Symbolics.get_variables-Tuple{DifferentialEquation}"},[s("span",{class:"jlbinding"},"Symbolics.get_variables")],-1)),i[7]||(i[7]=a()),t(n,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),i[8]||(i[8]=l('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(diff_eom</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Vector{Num}</span></span></code></pre></div><p>Return the dependent variables of <code>diff_eom</code>.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/afa425ceb096439a4ff601e5145532d557edf0e3/src/DifferentialEquation.jl#L25" target="_blank" rel="noreferrer">source</a></p>',3))]),s("details",g,[s("summary",null,[i[9]||(i[9]=s("a",{id:"HarmonicBalance.get_independent_variables-Tuple{DifferentialEquation}",href:"#HarmonicBalance.get_independent_variables-Tuple{DifferentialEquation}"},[s("span",{class:"jlbinding"},"HarmonicBalance.get_independent_variables")],-1)),i[10]||(i[10]=a()),t(n,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),i[11]||(i[11]=l(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_independent_variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    diff_eom</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DifferentialEquation</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Any</span></span></code></pre></div><p>Return the independent dependent variables of <code>diff_eom</code>.</p><p><a href="https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/afa425ceb096439a4ff601e5145532d557edf0e3/src/DifferentialEquation.jl#L43" target="_blank" rel="noreferrer">source</a></p>`,3))])])}const D=e(r,[["render",y]]);export{b as __pageData,D as default};
