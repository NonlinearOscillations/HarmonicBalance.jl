import{_ as l,c as t,a4 as n,j as i,a,o as e}from"./chunks/framework.7kS1RB1D.js";const p="/HarmonicBalance.jl/v0.10.5/assets/ohndoqr.DDG3oKdt.png",h="/HarmonicBalance.jl/v0.10.5/assets/qdsxmjc.Nsw0w518.png",k="/HarmonicBalance.jl/v0.10.5/assets/rgizkmj.CHo32oEM.png",T=JSON.parse('{"title":"Classifying solutions","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/classification.md","filePath":"tutorials/classification.md"}'),r={name:"tutorials/classification.md"},o={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},d={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.407ex",height:"1.027ex",role:"img",focusable:"false",viewBox:"0 -443 622 454","aria-hidden":"true"},g={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},c={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.027ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.319ex",height:"1.597ex",role:"img",focusable:"false",viewBox:"0 -694 583 706","aria-hidden":"true"},E={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},u={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.186ex"},xmlns:"http://www.w3.org/2000/svg",width:"8.359ex",height:"1.756ex",role:"img",focusable:"false",viewBox:"0 -694 3694.6 776","aria-hidden":"true"};function y(m,s,F,b,v,C){return e(),t("div",null,[s[11]||(s[11]=n(`<h1 id="classifying" tabindex="-1">Classifying solutions <a class="header-anchor" href="#classifying" aria-label="Permalink to &quot;Classifying solutions {#classifying}&quot;">​</a></h1><p>Given that you obtained some steady states for a parameter sweep of a specific model it can be useful to classify these solution. Let us consider a simple pametric oscillator</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> HarmonicBalance</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@variables</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω₀ γ λ α ω t </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(t)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">natural_equation </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, t), t) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> γ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, t) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (ω₀</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> λ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> cos</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> *</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ω </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> t)) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> α </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">^</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">diff_eq </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DifferentialEquation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(natural_equation, x)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add_harmonic!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(diff_eq, x, ω);</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">harmonic_eq </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> get_harmonic_equations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(diff_eq)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>A set of 2 harmonic equations</span></span>
<span class="line"><span>Variables: u1(T), v1(T)</span></span>
<span class="line"><span>Parameters: ω, α, γ, ω₀, λ</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Harmonic ansatz: </span></span>
<span class="line"><span>x(t) = u1(T)*cos(ωt) + v1(T)*sin(ωt)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Harmonic equations:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>-(1//2)*u1(T)*λ + (2//1)*Differential(T)(v1(T))*ω + Differential(T)(u1(T))*γ - u1(T)*(ω^2) + u1(T)*(ω₀^2) + v1(T)*γ*ω + (3//4)*(u1(T)^3)*α + (3//4)*u1(T)*(v1(T)^2)*α ~ 0</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Differential(T)(v1(T))*γ + (1//2)*v1(T)*λ - (2//1)*Differential(T)(u1(T))*ω - u1(T)*γ*ω - v1(T)*(ω^2) + v1(T)*(ω₀^2) + (3//4)*(u1(T)^2)*v1(T)*α + (3//4)*(v1(T)^3)*α ~ 0</span></span></code></pre></div>`,4)),i("p",null,[s[4]||(s[4]=a("We performe a 2d sweep in the driving frequency ")),i("mjx-container",o,[(e(),t("svg",d,s[0]||(s[0]=[i("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[i("g",{"data-mml-node":"math"},[i("g",{"data-mml-node":"mi"},[i("path",{"data-c":"1D714",d:"M495 384Q495 406 514 424T555 443Q574 443 589 425T604 364Q604 334 592 278T555 155T483 38T377 -11Q297 -11 267 66Q266 68 260 61Q201 -11 125 -11Q15 -11 15 139Q15 230 56 325T123 434Q135 441 147 436Q160 429 160 418Q160 406 140 379T94 306T62 208Q61 202 61 187Q61 124 85 100T143 76Q201 76 245 129L253 137V156Q258 297 317 297Q348 297 348 261Q348 243 338 213T318 158L308 135Q309 133 310 129T318 115T334 97T358 83T393 76Q456 76 501 148T546 274Q546 305 533 325T508 357T495 384Z",style:{"stroke-width":"3"}})])])],-1)]))),s[1]||(s[1]=i("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[i("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[i("mi",null,"ω")])],-1))]),s[5]||(s[5]=a(" and driving strength ")),i("mjx-container",g,[(e(),t("svg",c,s[2]||(s[2]=[i("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[i("g",{"data-mml-node":"math"},[i("g",{"data-mml-node":"mi"},[i("path",{"data-c":"1D706",d:"M166 673Q166 685 183 694H202Q292 691 316 644Q322 629 373 486T474 207T524 67Q531 47 537 34T546 15T551 6T555 2T556 -2T550 -11H482Q457 3 450 18T399 152L354 277L340 262Q327 246 293 207T236 141Q211 112 174 69Q123 9 111 -1T83 -12Q47 -12 47 20Q47 37 61 52T199 187Q229 216 266 252T321 306L338 322Q338 323 288 462T234 612Q214 657 183 657Q166 657 166 673Z",style:{"stroke-width":"3"}})])])],-1)]))),s[3]||(s[3]=i("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[i("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[i("mi",null,"λ")])],-1))]),s[6]||(s[6]=a(":"))]),s[12]||(s[12]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fixed </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (ω₀ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, γ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.002</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, α </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">varied </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (ω </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> range</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.99</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), λ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> range</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-6</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.03</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result_2D </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> get_steady_states</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(harmonic_eq, varied, fixed, threading</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>A steady state result for 10000 parameter points</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Solution branches:   5</span></span>
<span class="line"><span>   of which real:    5</span></span>
<span class="line"><span>   of which stable:  5</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Classes: stable, physical, Hopf, binary_labels</span></span></code></pre></div><p>By default the steady states of the system are classified by four different catogaries:</p><ul><li><p><code>physical</code>: Solutions that are physical, i.e., all variables are purely real.</p></li><li><p><code>stable</code>: Solutions that are stable, i.e., all eigenvalues of the Jacobian have negative real parts.</p></li><li><p><code>Hopf</code>: Solutions that are physical and have exactly two Jacobian eigenvalues with positive real parts, which are complex conjugates of each other. The class can help to identify regions where a limit cycle is present due to a <a href="https://en.wikipedia.org/wiki/Hopf_bifurcation" target="_blank" rel="noreferrer">Hopf bifurcation</a>. See also the tutorial on <a href="/HarmonicBalance.jl/v0.10.5/tutorials/limit_cycles#limit_cycles">limit cycles</a>.</p></li><li><p><code>binary_labels</code>: each region in the parameter sweep receives an identifier based on its permutation of stable branches. This allows to distinguish between different phases, which may have the same number of stable solutions.</p></li></ul><p>We can plot the number of stable solutions, giving the phase diagram</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><img src="`+p+'" alt=""></p>',7)),i("p",null,[s[9]||(s[9]=a("If we plot the a cut at ")),i("mjx-container",E,[(e(),t("svg",u,s[7]||(s[7]=[n('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="1D706" d="M166 673Q166 685 183 694H202Q292 691 316 644Q322 629 373 486T474 207T524 67Q531 47 537 34T546 15T551 6T555 2T556 -2T550 -11H482Q457 3 450 18T399 152L354 277L340 262Q327 246 293 207T236 141Q211 112 174 69Q123 9 111 -1T83 -12Q47 -12 47 20Q47 37 61 52T199 187Q229 216 266 252T321 306L338 322Q338 323 288 462T234 612Q214 657 183 657Q166 657 166 673Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(860.8,0)"><path data-c="3D" d="M56 347Q56 360 70 367H707Q722 359 722 347Q722 336 708 328L390 327H72Q56 332 56 347ZM56 153Q56 168 72 173H708Q722 163 722 153Q722 140 707 133H70Q56 140 56 153Z" style="stroke-width:3;"></path></g><g data-mml-node="mn" transform="translate(1916.6,0)"><path data-c="30" d="M96 585Q152 666 249 666Q297 666 345 640T423 548Q460 465 460 320Q460 165 417 83Q397 41 362 16T301 -15T250 -22Q224 -22 198 -16T137 16T82 83Q39 165 39 320Q39 494 96 585ZM321 597Q291 629 250 629Q208 629 178 597Q153 571 145 525T137 333Q137 175 145 125T181 46Q209 16 250 16Q290 16 318 46Q347 76 354 130T362 333Q362 478 354 524T321 597Z" style="stroke-width:3;"></path><path data-c="2E" d="M78 60Q78 84 95 102T138 120Q162 120 180 104T199 61Q199 36 182 18T139 0T96 17T78 60Z" transform="translate(500,0)" style="stroke-width:3;"></path><path data-c="30" d="M96 585Q152 666 249 666Q297 666 345 640T423 548Q460 465 460 320Q460 165 417 83Q397 41 362 16T301 -15T250 -22Q224 -22 198 -16T137 16T82 83Q39 165 39 320Q39 494 96 585ZM321 597Q291 629 250 629Q208 629 178 597Q153 571 145 525T137 333Q137 175 145 125T181 46Q209 16 250 16Q290 16 318 46Q347 76 354 130T362 333Q362 478 354 524T321 597Z" transform="translate(778,0)" style="stroke-width:3;"></path><path data-c="31" d="M213 578L200 573Q186 568 160 563T102 556H83V602H102Q149 604 189 617T245 641T273 663Q275 666 285 666Q294 666 302 660V361L303 61Q310 54 315 52T339 48T401 46H427V0H416Q395 3 257 3Q121 3 100 0H88V46H114Q136 46 152 46T177 47T193 50T201 52T207 57T213 61V578Z" transform="translate(1278,0)" style="stroke-width:3;"></path></g></g></g>',1)]))),s[8]||(s[8]=i("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[i("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[i("mi",null,"λ"),i("mo",null,"="),i("mn",null,"0.01")])],-1))]),s[10]||(s[10]=a(", we see that in the blue region only one stable solution exists with zero amplitude:"))]),s[13]||(s[13]=n(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, y</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;√(u1^2+v1^2)&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, cut</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">λ </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.01</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> display</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get_single_solution</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D; branch</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, index</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>OrderedCollections.OrderedDict{Num, ComplexF64} with 7 entries:</span></span>
<span class="line"><span>  u1 =&gt; -3.35208e-249-5.36333e-248im</span></span>
<span class="line"><span>  v1 =&gt; -7.59806e-248+1.45257e-248im</span></span>
<span class="line"><span>  ω  =&gt; 0.99+0.0im</span></span>
<span class="line"><span>  λ  =&gt; 1.0e-6+0.0im</span></span>
<span class="line"><span>  ω₀ =&gt; 1.0+0.0im</span></span>
<span class="line"><span>  γ  =&gt; 0.002+0.0im</span></span>
<span class="line"><span>  α  =&gt; 1.0+0.0im</span></span></code></pre></div><p>This solution becomes stable again outside the green lobe. Also called Mathieu lobe. Indeed, we can classify the zero amplitude solution by adding an extra catagory as a class:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">classify_solutions!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;sqrt(u1^2 + v1^2) &lt; 0.001&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;zero&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">result_2D</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>A steady state result for 10000 parameter points</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Solution branches:   5</span></span>
<span class="line"><span>   of which real:    5</span></span>
<span class="line"><span>   of which stable:  5</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Classes: zero, stable, physical, Hopf, binary_labels</span></span></code></pre></div><p>We can visualize the zero amplitude solution:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;zero&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p><img src="`+h+`" alt=""></p><p>This shows that inside the Mathieu lobe the zero amplitude solution becomes unstable due to the parametric drive being resonant with the oscillator.</p><p>We can also visualize the equi-amplitude curves of the solutions:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">classify_solutions!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;sqrt(u1^2 + v1^2) &gt; 0.12&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;large amplitude&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot_phase_diagram</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result_2D, class</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;large amplitude&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;stable&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p><img src="`+k+'" alt=""></p>',13))])}const x=l(r,[["render",y]]);export{T as __pageData,x as default};
