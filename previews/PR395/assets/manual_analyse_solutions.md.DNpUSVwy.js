import{_ as o,C as r,c as p,o as h,ai as n,j as a,G as t,a as e,w as l}from"./chunks/framework.05Bk7gqB.js";const C=JSON.parse('{"title":"Extracting Solutions","description":"","frontmatter":{},"headers":[],"relativePath":"manual/analyse_solutions.md","filePath":"manual/analyse_solutions.md"}'),c={name:"manual/analyse_solutions.md"},k={class:"jldocstring custom-block",open:""},d={class:"jldocstring custom-block",open:""},u={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""};function E(y,s,m,b,_,f){const i=r("Badge");return h(),p("div",null,[s[16]||(s[16]=n('<h1 id="Extracting-Solutions" tabindex="-1">Extracting Solutions <a class="header-anchor" href="#Extracting-Solutions" aria-label="Permalink to &quot;Extracting Solutions {#Extracting-Solutions}&quot;">​</a></h1><p>After computing the steady-states of the harmonic equations, you&#39;ll want to extract the solutions from the <a href="/HarmonicBalance.jl/previews/PR395/manual/solving_harmonics#HarmonicBalance.Result-manual-solving_harmonics"><code>HarmonicBalance.Result</code></a> struct.</p><h2 id="Basic-Solution-Extraction" tabindex="-1">Basic Solution Extraction <a class="header-anchor" href="#Basic-Solution-Extraction" aria-label="Permalink to &quot;Basic Solution Extraction {#Basic-Solution-Extraction}&quot;">​</a></h2><p>For plotting, you can extract the solutions using the <code>get_solutions</code> function, which parses a string into a symbolic expression, evaluates it for every steady state solution and filters the solutions by the requested class.</p>',4)),a("details",k,[a("summary",null,[s[0]||(s[0]=a("a",{id:"HarmonicBalance.get_solutions-manual-analyse_solutions",href:"#HarmonicBalance.get_solutions-manual-analyse_solutions"},[a("span",{class:"jlbinding"},"HarmonicBalance.get_solutions")],-1)),s[1]||(s[1]=e()),t(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[3]||(s[3]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_solutions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    branches</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">branch_count</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res), realify</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], not_class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    )</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_solutions</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; branches</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">branch_count</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res), class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], not_class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[])</span></span></code></pre></div><p>Extract solution vectors from a <code>Result</code> object based on specified filtering criteria given by the <code>class</code> keywords. The first method allows extracting a specific solution component by name <code>x</code>. The second method returns complete solution vectors.</p><p><strong>Keyword arguments</strong></p><ul><li><p><code>branches=1:branch_count(res)</code>: Range of branches to include in the output</p></li><li><p><code>realify=false</code>: Whether to convert complex solutions to real form</p></li><li><p><code>class=[&quot;physical&quot;, &quot;stable&quot;]</code>: Array of classification labels to include</p></li><li><p><code>not_class=[]</code>: Array of classification labels to exclude</p></li></ul><p><strong>Returns</strong></p><p>Filtered solution vectors matching the specified criteria</p>`,6)),t(i,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[2]||(s[2]=[a("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/82a1d77b7876da2f49a224a85fb6f146a98c8950/src/transform_solutions.jl#L194-L213",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),a("details",d,[a("summary",null,[s[4]||(s[4]=a("a",{id:"HarmonicBalance.get_single_solution-manual-analyse_solutions",href:"#HarmonicBalance.get_single_solution-manual-analyse_solutions"},[a("span",{class:"jlbinding"},"HarmonicBalance.get_single_solution")],-1)),s[5]||(s[5]=e()),t(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[7]||(s[7]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_single_solution</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">HarmonicBalance.Result{D, S, P, F}</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> F</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    branch,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    index</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Return an ordered dictionary specifying all variables and parameters of the solution in <code>result</code> on <code>branch</code> at the position <code>index</code>.</p>`,2)),t(i,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[6]||(s[6]=[a("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/82a1d77b7876da2f49a224a85fb6f146a98c8950/src/transform_solutions.jl#L1",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[17]||(s[17]=a("h2",{id:"attractors",tabindex:"-1"},[e("Attractors "),a("a",{class:"header-anchor",href:"#attractors","aria-label":'Permalink to "Attractors"'},"​")],-1)),a("details",u,[a("summary",null,[s[8]||(s[8]=a("a",{id:"HarmonicBalance.attractors-manual-analyse_solutions",href:"#HarmonicBalance.attractors-manual-analyse_solutions"},[a("span",{class:"jlbinding"},"HarmonicBalance.attractors")],-1)),s[9]||(s[9]=e()),t(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[11]||(s[11]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">attractors</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result{D}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, not_class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[]) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> D</span></span></code></pre></div><p>Extract attractors from a <a href="/HarmonicBalance.jl/previews/PR395/manual/solving_harmonics#HarmonicBalance.Result-manual-solving_harmonics"><code>Result</code></a> object. Returns an array of dictionaries, where each dictionary maps branch identifier to the attractor. The attractors are filtered by their corresponding class.</p><p><strong>Keyword arguments</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p><strong>Returns</strong></p><p><code>Array{Dict,D}</code>: Vector of dictionaries mapping branch indices to points satisfying the stability criteria at each parameter value</p>`,7)),t(i,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[10]||(s[10]=[a("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/82a1d77b7876da2f49a224a85fb6f146a98c8950/src/Result.jl#L122-L139",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),a("details",g,[a("summary",null,[s[12]||(s[12]=a("a",{id:"HarmonicBalance.phase_diagram-manual-analyse_solutions",href:"#HarmonicBalance.phase_diagram-manual-analyse_solutions"},[a("span",{class:"jlbinding"},"HarmonicBalance.phase_diagram")],-1)),s[13]||(s[13]=e()),t(i,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[15]||(s[15]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(res</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Result{D}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;physical&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, not_class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[]) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> {D}</span></span></code></pre></div><p>Calculate the phase diagram from a <code>Result</code> object by summing over the number of states at each swept parameters.</p><p><strong>Keyword arguments</strong></p><p>Class selection done by passing <code>String</code> or <code>Vector{String}</code> as kwarg:</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>class::String       :   only count solutions in this class (&quot;all&quot; --&gt; plot everything)</span></span>
<span class="line"><span>not_class::String   :   do not count solutions in this class</span></span></code></pre></div><p><strong>Returns</strong></p><ul><li>Array{Int64,D}: Sum of states after applying the specified class masks</li></ul>`,7)),t(i,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[14]||(s[14]=[a("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/82a1d77b7876da2f49a224a85fb6f146a98c8950/src/Result.jl#L92-L106",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})])])}const B=o(c,[["render",E]]);export{C as __pageData,B as default};
