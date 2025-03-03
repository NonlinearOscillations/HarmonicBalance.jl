import{_ as p,C as h,c as k,o,ai as n,j as i,G as e,a as t,w as l}from"./chunks/framework.DDVM-Dcu.js";const v=JSON.parse('{"title":"Time evolution","description":"","frontmatter":{},"headers":[],"relativePath":"manual/time_dependent.md","filePath":"manual/time_dependent.md"}'),r={name:"manual/time_dependent.md"},d={class:"jldocstring custom-block",open:""},E={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},y={class:"jldocstring custom-block",open:""};function u(m,s,b,F,C,f){const a=h("Badge");return o(),k("div",null,[s[24]||(s[24]=n("",3)),i("details",d,[i("summary",null,[s[0]||(s[0]=i("a",{id:"SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent",href:"#SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent"},[i("span",{class:"jlbinding"},"SciMLBase.ODEProblem")],-1)),s[1]||(s[1]=t()),e(a,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[3]||(s[3]=n("",2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[2]||(s[2]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/ext/TimeEvolution/ODEProblem.jl#L3-L9",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),i("details",E,[i("summary",null,[s[4]||(s[4]=i("a",{id:"HarmonicBalance.AdiabaticSweep-manual-time_dependent",href:"#HarmonicBalance.AdiabaticSweep-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.AdiabaticSweep")],-1)),s[5]||(s[5]=t()),e(a,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[7]||(s[7]=n("",12)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[6]||(s[6]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/src/types.jl#L9-L48",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[25]||(s[25]=i("p",null,[t("In addition, one can use the "),i("code",null,"steady_state_sweep"),t(" function from "),i("code",null,"SteadyStateDiffEqExt"),t(" to perform a parameter sweep over the steady states of a system. For this one has to load "),i("code",null,"SteadyStateDiffEq.jl"),t(".")],-1)),i("details",c,[i("summary",null,[s[8]||(s[8]=i("a",{id:"HarmonicBalance.steady_state_sweep-manual-time_dependent",href:"#HarmonicBalance.steady_state_sweep-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.steady_state_sweep")],-1)),s[9]||(s[9]=t()),e(a,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[12]||(s[12]=n("",2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[10]||(s[10]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/ext/SteadyStateDiffEqExt.jl#L12-L20",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1}),s[13]||(s[13]=n("",2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[11]||(s[11]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/ext/SteadyStateDiffEqExt.jl#L39-L47",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[26]||(s[26]=i("h2",{id:"plotting",tabindex:"-1"},[t("Plotting "),i("a",{class:"header-anchor",href:"#plotting","aria-label":'Permalink to "Plotting"'},"​")],-1)),i("details",g,[i("summary",null,[s[14]||(s[14]=i("a",{id:"RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent",href:"#RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent"},[i("span",{class:"jlbinding"},"RecipesBase.plot")],-1)),s[15]||(s[15]=t()),e(a,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),s[17]||(s[17]=n("",9)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[16]||(s[16]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/ext/PlotsExt/time_evolution.jl#L1-L19",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[27]||(s[27]=i("h2",{id:"miscellaneous",tabindex:"-1"},[t("Miscellaneous "),i("a",{class:"header-anchor",href:"#miscellaneous","aria-label":'Permalink to "Miscellaneous"'},"​")],-1)),s[28]||(s[28]=i("p",null,"Using a time-dependent simulation can verify solution stability in cases where the Jacobian is too expensive to compute.",-1)),i("details",y,[i("summary",null,[s[18]||(s[18]=i("a",{id:"HarmonicBalance.is_stable-manual-time_dependent",href:"#HarmonicBalance.is_stable-manual-time_dependent"},[i("span",{class:"jlbinding"},"HarmonicBalance.is_stable")],-1)),s[19]||(s[19]=t()),e(a,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),s[22]||(s[22]=n("",3)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[20]||(s[20]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/ext/TimeEvolution/ODEProblem.jl#L61",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1}),s[23]||(s[23]=n("",2)),e(a,{type:"info",class:"source-link",text:"source"},{default:l(()=>s[21]||(s[21]=[i("a",{href:"https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/96db5f27cc9d5221a4fe8fb2b363db02619524b6/src/classification.jl#L77",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})])])}const A=p(r,[["render",u]]);export{v as __pageData,A as default};
