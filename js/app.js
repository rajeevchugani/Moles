'use strict';
/* ---------- chemistry core ---------- */
const AR={H:1,C:12,N:14,O:16,Na:23,Mg:24,Al:27,P:31,S:32,Cl:35.5,K:39,Ca:40,Fe:56,Cu:64,Zn:65,Ag:108,Br:80};
const gcd=(a,b)=>b?gcd(b,a%b):a;
function parse(f){
  let i=0;
  const num=()=>{let s='';while(i<f.length&&/\d/.test(f[i]))s+=f[i++];return s?parseInt(s,10):1;};
  function grp(){
    const m={};
    while(i<f.length){
      const ch=f[i];
      if(ch==='('){i++;const sub=grp();i++;const n=num();for(const k in sub)m[k]=(m[k]||0)+sub[k]*n;}
      else if(ch===')')return m;
      else if(/[A-Z]/.test(ch)){let el=ch;i++;while(i<f.length&&/[a-z]/.test(f[i]))el+=f[i++];const n=num();m[el]=(m[el]||0)+n;}
      else i++;
    }
    return m;
  }
  return grp();
}
const mr=f=>{const p=parse(f);let s=0;for(const e in p)s+=(AR[e]||0)*p[e];return s;};
const fx=f=>String(f).replace(/([A-Za-z\)])(\d+)/g,'$1<sub>$2</sub>');
const sub=f=>String(f).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[d]);
const chg=(c,s)=>(c===1?'':c)+s;
const fmt=x=>{if(x===0)return '0';const v=parseFloat(Number(x).toPrecision(4));return String(v);};
const norm=s=>String(s).replace(/[₀₁₂₃₄₅₆₇₈₉]/g,d=>'₀₁₂₃₄₅₆₇₈₉'.indexOf(d)).replace(/[\s_{}^]/g,'').replace(/\[/g,'(').replace(/\]/g,')');

const CAT={
  Na:{f:'Na',c:1,n:'sodium'},K:{f:'K',c:1,n:'potassium'},Mg:{f:'Mg',c:2,n:'magnesium'},Ca:{f:'Ca',c:2,n:'calcium'},
  Al:{f:'Al',c:3,n:'aluminium'},Zn:{f:'Zn',c:2,n:'zinc'},Fe2:{f:'Fe',c:2,n:'iron(II)'},Fe3:{f:'Fe',c:3,n:'iron(III)'},
  Cu:{f:'Cu',c:2,n:'copper(II)'},Ag:{f:'Ag',c:1,n:'silver'},NH4:{f:'NH4',c:1,n:'ammonium',poly:1}
};
const AN={
  O:{f:'O',c:2,n:'oxide'},Cl:{f:'Cl',c:1,n:'chloride'},S:{f:'S',c:2,n:'sulfide'},SO4:{f:'SO4',c:2,n:'sulfate',poly:1},
  NO3:{f:'NO3',c:1,n:'nitrate',poly:1},PO4:{f:'PO4',c:3,n:'phosphate',poly:1},CO3:{f:'CO3',c:2,n:'carbonate',poly:1},
  OH:{f:'OH',c:1,n:'hydroxide',poly:1},Ac:{f:'CH3COO',c:1,n:'ethanoate',poly:1}
};
function salt(ck,ak){
  const c=CAT[ck],a=AN[ak],g=gcd(c.c,a.c),nc=a.c/g,na=c.c/g;
  const part=(x,n)=>n===1?x.f:(x.poly?`(${x.f})${n}`:x.f+n);
  const f=(ak==='Ac'&&nc===1&&na===1)?'CH3COO'+c.f:part(c,nc)+part(a,na);
  return {f,name:`${c.n} ${a.n}`,ion:{c,a}};
}
const soluble=(ck,ak)=>['Na','K','NH4'].includes(ck)||['Cl','NO3','Ac'].includes(ak)||(ak==='SO4'&&ck!=='Ca');
const MN={Na:'sodium',K:'potassium',Mg:'magnesium',Ca:'calcium',Al:'aluminium',Zn:'zinc',Fe:'iron',Cu:'copper',Ag:'silver'};
const SP=(f,name,form,x)=>Object.assign({f,name,form},x||{});
const wname=s=>s.form==='aq'&&!/acid|solution/.test(s.name)?s.name+' solution':s.name;

/* ---------- reaction library ---------- */
const CATS={
  synthesis:{label:'Synthesis',pattern:'metal + non-metal → ionic compound',info:'Two elements join to make one product. A metal and a non-metal make an ionic compound, named from the metal then the non-metal with an -ide ending.'},
  combustion:{label:'Combustion',pattern:'fuel + oxygen → oxide(s)',info:'A fuel burns in oxygen and releases energy. With plenty of oxygen carbon becomes carbon dioxide and hydrogen becomes water. With limited oxygen carbon monoxide forms.'},
  decomposition:{label:'Decomposition',pattern:'one compound → simpler substances',info:'One compound is split by heat (thermal decomposition), electricity or a catalyst into two or more simpler substances.'},
  displacement:{label:'Displacement',pattern:'metal + compound → new metal (or H₂) + new compound',info:'A more reactive metal pushes a less reactive metal out of its compound. Metals above hydrogen in the reactivity series also displace hydrogen from dilute acids.'},
  neutralisation:{label:'Acid reactions',pattern:'acid + base → salt + water',info:'Acids react with metals (salt + hydrogen), metal oxides and hydroxides (salt + water) and metal carbonates (salt + water + carbon dioxide). The salt name comes from the acid: chloride, sulfate, nitrate, phosphate, carbonate, ethanoate.'}
};
const MET_SYN=['Na','K','Mg','Ca','Al','Zn','Fe','Cu'];
const FUELS={H2:['hydrogen','gas'],C:['carbon','solid'],S:['sulfur','solid'],CH4:['methane','gas'],C2H6:['ethane','gas'],C3H8:['propane','gas'],C4H10:['butane','gas'],C2H4:['ethene','gas'],CH3OH:['methanol','liquid'],C2H5OH:['ethanol','liquid']};
const DEC={
  CaCO3:{r:['CaCO3','calcium carbonate','solid'],p:[['CaO','calcium oxide','solid'],['CO2','carbon dioxide','gas']],cond:'heat'},
  CuCO3:{r:['CuCO3','copper(II) carbonate','solid'],p:[['CuO','copper(II) oxide','solid'],['CO2','carbon dioxide','gas']],cond:'heat'},
  ZnCO3:{r:['ZnCO3','zinc carbonate','solid'],p:[['ZnO','zinc oxide','solid'],['CO2','carbon dioxide','gas']],cond:'heat'},
  MgCO3:{r:['MgCO3','magnesium carbonate','solid'],p:[['MgO','magnesium oxide','solid'],['CO2','carbon dioxide','gas']],cond:'heat'},
  'Cu(OH)2':{r:['Cu(OH)2','copper(II) hydroxide','solid'],p:[['CuO','copper(II) oxide','solid'],['H2O','water','liquid']],cond:'heat'},
  'Cu(NO3)2':{r:['Cu(NO3)2','copper(II) nitrate','solid'],p:[['CuO','copper(II) oxide','solid'],['NO2','nitrogen dioxide','gas'],['O2','oxygen','gas']],cond:'heat'},
  'Zn(NO3)2':{r:['Zn(NO3)2','zinc nitrate','solid'],p:[['ZnO','zinc oxide','solid'],['NO2','nitrogen dioxide','gas'],['O2','oxygen','gas']],cond:'heat'},
  NaNO3:{r:['NaNO3','sodium nitrate','solid'],p:[['NaNO2','sodium nitrite','solid'],['O2','oxygen','gas']],cond:'heat'},
  KNO3:{r:['KNO3','potassium nitrate','solid'],p:[['KNO2','potassium nitrite','solid'],['O2','oxygen','gas']],cond:'heat'},
  NaHCO3:{r:['NaHCO3','sodium hydrogencarbonate','solid'],p:[['Na2CO3','sodium carbonate','solid'],['H2O','water','liquid'],['CO2','carbon dioxide','gas']],cond:'heat'},
  H2O2:{r:['H2O2','hydrogen peroxide','aq'],p:[['H2O','water','liquid'],['O2','oxygen','gas']],cond:'MnO₂ catalyst'},
  H2O:{r:['H2O','water','liquid'],p:[['H2','hydrogen','gas'],['O2','oxygen','gas']],cond:'electrolysis'}
};
const RANK={Mg:6,Al:5,Zn:4,Fe:3,H:2,Cu:1,Ag:0};
const DMET={Mg:'Mg',Al:'Al',Zn:'Zn',Fe:'Fe2',Cu:'Cu'};
const DPART={
  CuSO4:{m:'Cu',a:'SO4'},FeSO4:{m:'Fe',a:'SO4'},ZnSO4:{m:'Zn',a:'SO4'},MgSO4:{m:'Mg',a:'SO4'},
  AgNO3:{m:'Ag',a:'NO3'},'Cu(NO3)2':{m:'Cu',a:'NO3'},'Zn(NO3)2':{m:'Zn',a:'NO3'},'Mg(NO3)2':{m:'Mg',a:'NO3'},
  HCl:{acid:1,a:'Cl'},H2SO4:{acid:1,a:'SO4'}
};
const ACIDS={HCl:['hydrochloric acid','Cl'],H2SO4:['sulfuric acid','SO4'],H3PO4:['phosphoric acid','PO4'],HNO3:['nitric acid','NO3'],H2CO3:['carbonic acid','CO3'],CH3COOH:['ethanoic acid','Ac']};
const BASES={
  NaOH:['sodium hydroxide','hyd','Na','aq'],KOH:['potassium hydroxide','hyd','K','aq'],'Ca(OH)2':['calcium hydroxide','hyd','Ca','aq'],'Mg(OH)2':['magnesium hydroxide','hyd','Mg','solid'],
  MgO:['magnesium oxide','oxide','Mg','solid'],CuO:['copper(II) oxide','oxide','Cu','solid'],ZnO:['zinc oxide','oxide','Zn','solid'],CaO:['calcium oxide','oxide','Ca','solid'],
  Na2CO3:['sodium carbonate','carb','Na','aq'],K2CO3:['potassium carbonate','carb','K','aq'],CaCO3:['calcium carbonate','carb','Ca','solid'],
  CuCO3:['copper(II) carbonate','carb','Cu','solid'],ZnCO3:['zinc carbonate','carb','Zn','solid'],MgCO3:['magnesium carbonate','carb','Mg','solid'],
  NH3:['ammonia','nh3','NH4','aq']
};
const NTYPE={metal:{label:'Metal + acid',pat:'metal + acid → salt + hydrogen'},oxide:{label:'Metal oxide + acid',pat:'metal oxide + acid → salt + water'},hyd:{label:'Metal hydroxide + acid',pat:'metal hydroxide + acid → salt + water'},carb:{label:'Metal carbonate + acid',pat:'metal carbonate + acid → salt + water + carbon dioxide'},nh3:{label:'Ammonia + acid',pat:'ammonia + acid → ammonium salt'}};
const NBASE={metal:['Mg','Al','Zn','Fe','Cu'],oxide:['MgO','CuO','ZnO','CaO'],hyd:['NaOH','KOH','Ca(OH)2','Mg(OH)2'],carb:['Na2CO3','K2CO3','CaCO3','CuCO3','ZnCO3','MgCO3'],nh3:['NH3']};
const nbName=k=>MN[k]?`${MN[k]} (${k})`:`${BASES[k][0]} (${sub(k)})`;
function optsFor(cat,p){
  if(cat!=='neutralisation')return OPTS[cat];
  const t=NBASE[p.type]?p.type:'hyd';
  return [{k:'type',label:'Type of reaction',list:Object.keys(NTYPE).map(k=>[k,NTYPE[k].label])},
          {k:'acid',label:'Acid',list:OPTS.neutralisation[0].list},
          {k:'base',label:{metal:'Metal',oxide:'Metal oxide',hyd:'Metal hydroxide',carb:'Metal carbonate',nh3:'Base'}[t],list:NBASE[t].map(k=>[k,nbName(k)])}];
}
const OPTS={
  synthesis:[{k:'metal',label:'Metal',list:MET_SYN.map(m=>[m,MN[m]+` (${m})`])},{k:'nonmetal',label:'Non-metal',list:[['O2','oxygen (O₂)'],['Cl2','chlorine (Cl₂)'],['S','sulfur (S)']]}],
  combustion:[{k:'fuel',label:'Fuel',list:Object.keys(FUELS).map(f=>[f,`${FUELS[f][0]} (${sub(f)})`])},{k:'oxy',label:'Oxygen supply',list:[['complete','Plenty of oxygen'],['incomplete','Limited oxygen']]}],
  decomposition:[{k:'dec',label:'Compound to split',list:Object.keys(DEC).map(k=>[k,`${DEC[k].r[1]} (${sub(k)})`])}],
  displacement:[{k:'metal',label:'Metal added',list:Object.keys(DMET).map(m=>[m,`${MN[m]} (${m})`])},{k:'partner',label:'Solution it goes into',list:Object.keys(DPART).map(k=>[k,(DPART[k].acid?ACIDS[k][0]:salt(DPART[k].m==='Fe'?'Fe2':DPART[k].m,DPART[k].a).name+' solution')+` (${sub(k)})`])}],
  neutralisation:[{k:'acid',label:'Acid',list:Object.keys(ACIDS).map(k=>[k,`${ACIDS[k][0]} (${sub(k)})`])},{k:'base',label:'Base, alkali or carbonate',list:Object.keys(BASES).map(k=>[k,`${BASES[k][0]} (${sub(k)})`])}]
};
const DEF_PICK={synthesis:{metal:'Mg',nonmetal:'O2'},combustion:{fuel:'CH4',oxy:'complete'},decomposition:{dec:'CaCO3'},displacement:{metal:'Mg',partner:'CuSO4'},neutralisation:{type:'hyd',acid:'HCl',base:'NaOH'}};

function build(cat,p){
  const out={cat,R:[],P:[],cond:'',valid:true,msg:''};
  if(cat==='synthesis'){
    const nm=p.metal===undefined?'':p.nonmetal;
    const an={O2:'O',Cl2:'Cl',S:'S'}[nm];
    const ck=p.metal==='Fe'?(nm==='S'?'Fe2':'Fe3'):p.metal;
    const s=salt(ck,an);
    out.R=[SP(p.metal,MN[p.metal],'solid'),nm==='O2'?SP('O2','oxygen','gas'):nm==='Cl2'?SP('Cl2','chlorine','gas'):SP('S','sulfur','solid')];
    out.P=[SP(s.f,s.name,'solid',{ion:s.ion})];out.cond='heat';
  }else if(cat==='combustion'){
    const [n,form]=FUELS[p.fuel],el=parse(p.fuel),inc=p.oxy==='incomplete'&&el.C;
    out.R=[SP(p.fuel,n,form),SP('O2','oxygen','gas')];
    if(el.C)out.P.push(inc?SP('CO','carbon monoxide','gas'):SP('CO2','carbon dioxide','gas'));
    if(el.H)out.P.push(SP('H2O','water','liquid'));
    if(el.S)out.P.push(SP('SO2','sulfur dioxide','gas'));
    out.cond='burn';
  }else if(cat==='decomposition'){
    const d=DEC[p.dec];out.R=[SP(...d.r)];out.P=d.p.map(x=>SP(...x));out.cond=d.cond;
  }else if(cat==='displacement'){
    const part=DPART[p.partner],mk=DMET[p.metal];
    out.R=[SP(p.metal,MN[p.metal],'solid')];
    const rp=part.acid?RANK.H:RANK[part.m];
    if(part.acid)out.R.push(SP(p.partner,ACIDS[p.partner][0],'aq'));
    else{const s=salt(part.m==='Fe'?'Fe2':part.m,part.a);out.R.push(SP(s.f,s.name,'aq',{ion:s.ion}));}
    if(RANK[p.metal]<=rp){out.valid=false;out.msg=part.acid?`${MN[p.metal]} is below hydrogen in the reactivity series, so it cannot displace hydrogen from the acid. No reaction.`:`${MN[p.metal]} is not more reactive than ${MN[part.m]}, so it cannot displace it. No reaction.`;}
    else{
      const s=salt(mk,part.a);
      out.P=[SP(s.f,s.name,'aq',{ion:s.ion}),part.acid?SP('H2','hydrogen','gas'):SP(part.m,MN[part.m],'solid')];
    }
  }else if(p.type==='metal'){
    const [an,ak]=ACIDS[p.acid],mk=DMET[p.base];
    out.type='metal';
    out.R=[SP(p.acid,an,'aq'),SP(p.base,MN[p.base],'solid')];
    if(p.acid==='HNO3'){out.valid=false;out.msg='Nitric acid is an oxidising acid. With metals it gives nitrogen oxides instead of hydrogen, so it is not used here. Choose another acid.';}
    else if(p.acid==='H2CO3'){out.valid=false;out.msg='Carbonic acid is too weak to react with a metal at this level. Choose another acid.';}
    else if(RANK[p.base]<=RANK.H){out.valid=false;out.msg=`${MN[p.base]} is below hydrogen in the reactivity series, so it does not react with dilute acids. No reaction.`;}
    else{const s=salt(mk,ak);out.P=[SP(s.f,s.name,soluble(mk==='Fe2'?'Fe2':mk,ak)?'aq':'solid',{ion:s.ion}),SP('H2','hydrogen','gas')];}
  }else{
    out.type=p.type;
    const [an,ak]=ACIDS[p.acid],[bn,bt,bc,bf]=BASES[p.base];
    out.R=[SP(p.acid,an,'aq'),SP(p.base,bn,bf)];
    if(p.acid==='H2CO3'&&bt==='carb'){out.valid=false;out.msg='Carbonic acid is itself a carbonate source and is too weak to react with a carbonate in this way. Choose a different acid or base.';}
    else{
      const s=salt(bc,ak),sol=soluble(bc,ak);
      out.P=[SP(s.f,s.name,sol?'aq':'solid',{ion:s.ion})];
      if(bt!=='nh3')out.P.push(SP('H2O','water','liquid'));
      if(bt==='carb')out.P.push(SP('CO2','carbon dioxide','gas'));
    }
  }
  out.R.forEach(s=>s.name=wname(s));out.P.forEach(s=>s.name=wname(s));
  if(out.valid){
    const co=balance(out.R.map(s=>s.f),out.P.map(s=>s.f));
    if(!co){out.valid=false;out.msg='This combination cannot be balanced as shown. Try another.';}else out.coefs=co;
  }
  return out;
}

/* balancing: solve atom conservation as a null-space problem */
function balance(rs,ps){
  const all=[...rs,...ps],parsed=all.map(parse),els=[...new Set(parsed.flatMap(o=>Object.keys(o)))];
  const M=els.map(e=>parsed.map((o,j)=>(o[e]||0)*(j<rs.length?1:-1)));
  const rows=M.length,cols=all.length,piv=[];let r=0;
  for(let c=0;c<cols&&r<rows;c++){
    let b=r;for(let i=r+1;i<rows;i++)if(Math.abs(M[i][c])>Math.abs(M[b][c]))b=i;
    if(Math.abs(M[b][c])<1e-9)continue;
    [M[r],M[b]]=[M[b],M[r]];
    const d=M[r][c];for(let j=0;j<cols;j++)M[r][j]/=d;
    for(let i=0;i<rows;i++)if(i!==r){const k=M[i][c];if(k)for(let j=0;j<cols;j++)M[i][j]-=k*M[r][j];}
    piv.push(c);r++;
  }
  const free=[];for(let c=0;c<cols;c++)if(!piv.includes(c))free.push(c);
  if(free.length!==1)return null;
  const v=new Array(cols).fill(0);v[free[0]]=1;
  piv.forEach((c,i)=>{v[c]=-M[i][free[0]];});
  for(let k=1;k<=120;k++){
    if(v.every(x=>Math.abs(x*k-Math.round(x*k))<1e-6&&x*k>0.5))return v.map(x=>Math.round(x*k));
  }
  return null;
}

/* ---------- amounts and calculations ---------- */
const nOf=(s,a)=>s.form==='aq'?a.c*a.V/1000:s.form==='gas'?a.V/24:a.m/mr(s.f);
const unitOf=s=>s.form==='gas'?'dm³':'g';
const amtOf=(s,n)=>s.form==='gas'?n*24:n*mr(s.f);
function calc(R,A){
  const nr=R.R.length,co=R.coefs;
  const n=R.R.map((s,i)=>nOf(s,A[i]));
  if(n.some(x=>!(x>0)))return {bad:true};
  const ratio=n.map((x,i)=>x/co[i]),xi=Math.min(...ratio);
  const lim=ratio.map(r=>r-xi<=xi*0.005);
  const used=co.slice(0,nr).map(c=>c*xi),left=n.map((x,i)=>Math.max(0,x-used[i]));
  const nP=R.P.map((s,j)=>co[nr+j]*xi);
  const Vt=R.R.reduce((t,s,i)=>t+(s.form==='aq'?A[i].V:0),0);
  return {n,ratio,xi,lim,used,left,leftAmt:R.R.map((s,i)=>amtOf(s,left[i])),nP,pAmt:R.P.map((s,j)=>amtOf(s,nP[j])),Vt,
    pConc:R.P.map((s,j)=>s.form==='aq'&&Vt>0?nP[j]/(Vt/1000):null),allLim:lim.every(Boolean)};
}
const rnd=a=>a[Math.floor(Math.random()*a.length)];
function defaults(R){
  const co=R.coefs,nr=R.R.length;
  for(let tries=0;tries<30;tries++){
    const base=rnd([0.02,0.025,0.04,0.05,0.075,0.1]),A=[];
    let ok=true;
    R.R.forEach((s,i)=>{
      let f=i===0?1:rnd([0.5,0.7,1.3,1.6,2]);
      const n=base*co[i]/co[0]*f;
      if(s.form==='solid'||s.form==='liquid'){const m=Math.round(n*mr(s.f)*10)/10;if(m<0.4||m>60)ok=false;A.push({m});}
      else if(s.form==='gas'){const V=Math.round(n*24*10)/10;if(V<0.2||V>15)ok=false;A.push({V});}
      else{const c=rnd([0.1,0.2,0.5,1,2]);const V=Math.max(10,Math.round(n/c*1000/5)*5);if(V>250)ok=false;A.push({c,V});}
    });
    if(ok)return A;
  }
  return R.R.map(s=>s.form==='aq'?{c:1,V:50}:s.form==='gas'?{V:2.4}:{m:5});
}

/* ---------- state ---------- */
const MODES=[['worked','Worked','All the working is shown. Step through it and read the reasons.'],
 ['guided','Guided','Hints stay on and every answer is marked as you type. Each calculation step unlocks when the last one is right.'],
 ['standard','Standard','More blanks. Hints on request and you press Check when ready.'],
 ['challenge','Challenge','Almost everything is a blank, including coefficients and which formula to use. No hints until you press Check.']];
const S={mode:'guided',tab:'lab',step:1,cat:'synthesis',pick:JSON.parse(JSON.stringify(DEF_PICK)),R:null,A:[],vals:{},checked:{},shown:{},fields:{},coefs:[],hintOpen:{},okf:{},prac:'mix',tit:null,run:null,reveal:0};
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;');
const allSp=()=>S.R.R.concat(S.R.P);
const supd=s=>String(s).replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);

function setReaction(){
  S.R=build(S.cat,S.pick[S.cat]);
  S.vals={};S.checked={};S.shown={};S.hintOpen={};S.okf={};S.run=null;S.tit=null;S.prac='mix';S.reveal=0;
  S.coefs=S.R.valid?allSp().map(()=>1):[];
  if(S.R.valid&&S.mode==='worked')S.coefs=S.R.coefs.slice();
  S.A=S.R.valid?defaults(S.R):[];
  renderAll();
}
function setMode(m){
  S.mode=m;S.vals={};S.checked={};S.shown={};S.hintOpen={};S.okf={};S.reveal=0;
  if(S.R&&S.R.valid)S.coefs=m==='worked'?S.R.coefs.slice():allSp().map(()=>1);
  renderAll();
}

/* ---------- fields: blanks that can be checked ---------- */
function inp(id,ans,o){
  o=o||{};S.fields[id]={ans,kind:o.kind||'num',tol:o.tol};
  const v=S.shown[id]?ansStr(id):(S.vals[id]!==undefined?S.vals[id]:'');
  return `<input class="blank${o.wide?' wide':''}" data-f="${id}" type="text" inputmode="${o.kind==='txt'?'text':'decimal'}" autocomplete="off" spellcheck="false" value="${esc(v)}" aria-label="${esc(o.label||id)}"${o.dis?' disabled':''}>`;
}
function selF(id,ans,options,o){
  o=o||{};S.fields[id]={ans,kind:'sel'};
  const v=S.shown[id]?ans:(S.vals[id]||'');
  return `<select class="blank" data-f="${id}" aria-label="${esc(o.label||id)}"${o.dis?' disabled':''}>${options.map(([k,l])=>`<option value="${k}"${k===v?' selected':''}>${l}</option>`).join('')}</select>`;
}
const ansStr=id=>{const f=S.fields[id];return f.kind==='num'?fmt(f.ans):f.ans;};
function sameF(v,ans){
  const a=norm(v);if(!/^[A-Za-z0-9()]+$/.test(a))return false;
  if(a===norm(ans))return true;
  const p=parse(a),q=parse(ans),ks=Object.keys(q);
  return Object.keys(p).length===ks.length&&ks.every(k=>p[k]===q[k]);
}
function isRight(id){
  const f=S.fields[id],v=S.vals[id];
  if(!f||v===undefined||v==='')return null;
  if(f.kind==='txt')return sameF(v,f.ans);
  if(f.kind==='sel')return v===f.ans;
  const x=parseFloat(String(v).replace(',','.'));
  if(isNaN(x))return false;
  const a=f.ans;
  if(Math.abs(a)<1e-9)return Math.abs(x)<0.0005;
  return Math.abs(x-a)<=Math.abs(a)*(f.tol||0.02)+1e-9;
}
function markOf(id){
  if(S.shown[id])return 'shown';
  const r=isRight(id);if(r===null)return '';
  return (S.mode==='guided'||S.checked[id])?(r?'ok':'bad'):'';
}
function paint(){
  Object.keys(S.fields).forEach(id=>{if(/^f\d+$/.test(id)&&isRight(id)===true&&(S.mode==='guided'||S.checked[id]))S.okf[id]=true;});
  document.querySelectorAll('[data-f]').forEach(el=>{
    const m=markOf(el.dataset.f);el.classList.remove('ok','bad','shown');if(m)el.classList.add(m);
  });
  if(typeof updateLocks==='function')updateLocks();
  renderEq();renderScore();
}
function renderScore(){
  const el=$('score');if(!el)return;
  const ids=Object.keys(S.fields);let ok=0;
  ids.forEach(id=>{if(isRight(id)===true&&!S.shown[id])ok++;});
  el.textContent=ids.length?`${ok} of ${ids.length} blanks correct`:'';
}
const fok=i=>S.mode==='worked'||!!S.okf['f'+i]||!!S.shown['f'+i];
function balInfo(){
  const sp=allSp(),nr=S.R.R.length,els=[...new Set(sp.flatMap(s=>Object.keys(parse(s.f))))];
  return els.map(e=>{let l=0,r=0;sp.forEach((s,i)=>{const c=(parse(s.f)[e]||0)*S.coefs[i];if(i<nr)l+=c;else r+=c;});return {e,l,r,ok:l===r};});
}
const isBal=()=>S.R&&S.R.valid&&S.coefs.length&&balInfo().every(x=>x.ok);
const simplest=()=>S.coefs.reduce((a,b)=>gcd(a,b))===1;

/* ---------- header, rail, equation bar ---------- */
function renderModes(){
  $('modeSeg').innerHTML=MODES.map(([k,l])=>`<button data-act="mode" data-m="${k}" aria-pressed="${S.mode===k}">${l}</button>`).join('');
  $('modeNote').textContent=MODES.find(m=>m[0]===S.mode)[2];
  $('tabs').innerHTML=[['lab','Reaction lab'],['tool','Mole practice'],['bal','Balancing practice'],['yield','% yield practice']].map(([k,l])=>`<button data-act="tab" data-t="${k}" aria-pressed="${S.tab===k}">${l}</button>`).join('');
  $('labView').hidden=S.tab!=='lab';$('toolView').hidden=S.tab!=='tool';$('balView').hidden=S.tab!=='bal';$('yView').hidden=S.tab!=='yield';
}
const STEPS=['Reaction','Formulae','Balance','Experiment','Calculate'];
function renderRail(){
  const R=S.R,v=R&&R.valid;
  const done=[v,v&&allSp().every((s,i)=>fok(i)),v&&isBal()&&(S.mode==='worked'||true),!!S.run||(S.tit&&S.tit.done),false];
  $('rail').innerHTML=STEPS.map((l,i)=>`<li><button data-act="go" data-s="${i+1}" ${S.step===i+1?'aria-current="step"':''} class="${done[i]&&S.step!==i+1?'done':''}"><b>${done[i]&&S.step!==i+1?'✓':i+1}</b>${l}</button></li>`).join('');
}
function renderEq(){
  const R=S.R;if(!R)return;
  const sp=allSp(),nr=R.R.length;
  const word=l=>l.map(s=>`<span class="sp">${s.name}</span>`).join('<span class="plus">+</span>');
  const arrow=`<span class="arr">${R.cond||'&nbsp;'}<i>→</i></span>`;
  let h=`<div class="lbl">Word equation</div><div class="eq">${word(R.R)}${R.valid?arrow+word(R.P):'<span class="arr">&nbsp;<i>↛</i></span><span class="sp">no reaction</span>'}</div>`;
  if(R.valid){
    const showCo=S.mode==='worked'||(S.step>=3&&isBal());
    const co=S.mode==='worked'?R.coefs:S.coefs;
    const st={solid:'s',liquid:'l',gas:'g',aq:'aq'};
    const chem=(s,i)=>{
      const known=fok(i)||S.step>=3&&S.mode==='worked';
      const c=showCo&&co[i]>1?co[i]+' ':'';
      return known?`<span class="sp">${c}${fx(s.f)}<span class="st">(${st[s.form]})</span></span>`:`<span class="sp unk" aria-label="not yet written">?</span>`;
    };
    h+=`<div class="lbl" style="margin-top:4px">Chemical equation${showCo?' (balanced)':''}</div><div class="eq chem">${R.R.map((s,i)=>chem(s,i)).join('<span class="plus">+</span>')}${arrow}${R.P.map((s,j)=>chem(s,nr+j)).join('<span class="plus">+</span>')}</div>`;
  }
  $('eqbar').innerHTML=h;
}
function navHTML(next){
  const nx=S.step<5;
  return `<div class="nav"><span>${S.step>1?`<button class="btn" data-act="go" data-s="${S.step-1}">← Back</button>`:''}</span><span class="score" id="score"></span><span>${nx?`<button class="btn pri" data-act="go" data-s="${S.step+1}"${next===false?' disabled':''}>Next: ${STEPS[S.step]} →</button>`:''}</span></div>`;
}
function renderStage(){
  S.fields={};
  [null,stage1,stage2,stage3,stage4,stage5][S.step]();
  renderScore();
}
function renderAll(){
  renderModes();
  if(S.tab==='lab'){renderRail();renderEq();renderStage();}else if(S.tab==='tool')renderTool();else if(S.tab==='bal')renderBal();else renderYield();
}
function go(n){S.step=n;renderRail();renderEq();renderStage();$('rail').scrollIntoView({block:'start',behavior:'smooth'});}

/* ---------- step 1: choose ---------- */
function stage1(){
  const R=S.R,c=CATS[S.cat],p=S.pick[S.cat];
  const noC=S.cat==='combustion'&&!parse(p.fuel).C;
  const sels=optsFor(S.cat,p).map(o=>`<label class="f">${o.label}<select data-pick="${o.k}"${o.k==='oxy'&&noC?' disabled':''}>${o.list.map(([k,l])=>`<option value="${k}"${p[o.k]===k?' selected':''}>${l}</option>`).join('')}</select></label>`).join('');
  const rs=S.cat==='displacement'?`<details><summary>Reactivity series (most to least reactive)</summary><p style="margin-top:8px" class="given">K &gt; Na &gt; Ca &gt; Mg &gt; Al &gt; Zn &gt; Fe &gt; (H) &gt; Cu &gt; Ag</p></details>`:'';
  $('stage').innerHTML=`<header><h2>Step 1 · Choose the reaction</h2><p>Pick a type of reaction, then choose your own reactants. The word equation above updates straight away.</p></header>
  <div class="chips" role="group" aria-label="Reaction type">${Object.keys(CATS).map(k=>`<button class="chip" data-act="cat" data-c="${k}" aria-pressed="${S.cat===k}">${CATS[k].label}</button>`).join('')}</div>
  <div class="callout"><b>${S.cat==='neutralisation'&&NTYPE[p.type]?NTYPE[p.type].pat:c.pattern}</b><br>${c.info}</div>
  <div class="pick">${sels}</div>
  <div class="chips"><button class="btn" data-act="random">Random reaction</button></div>
  ${rs}
  ${R.valid?`<div class="callout good">Reaction ready. Next, work out the chemical formula for every substance in the word equation.</div>`:`<div class="callout bad">${R.msg}</div>`}
  ${navHTML(R.valid)}`;
}

/* ---------- step 2: formulae ---------- */
function hintF(s){
  if(s.ion){const c=s.ion.c,a=s.ion.a;return `<b>Positive ion:</b> ${sub(c.f)}${supd(c.c===1?'':c.c)}⁺ (${c.n})<br><b>Negative ion:</b> ${sub(a.f)}${supd(a.c===1?'':a.c)}⁻ (${a.n})<br>Balance the charges so they add to zero.`;}
  if(['H2','O2','Cl2','N2'].includes(s.f))return 'A diatomic element: its molecules are two atoms joined together.';
  if(Object.keys(parse(s.f)).length===1&&!/\d/.test(s.f))return 'An element on its own: just its symbol.';
  const n=Object.values(parse(s.f)).reduce((a,b)=>a+b,0);
  return `Starts with “${s.f[0]}” and has ${n} atoms in each molecule or formula unit.`;
}
function stage2(){
  const R=S.R;
  if(!R.valid){$('stage').innerHTML=`<header><h2>Step 2 · Formulae</h2></header><div class="callout bad">${R.msg}</div>${navHTML(false)}`;return;}
  const sp=allSp(),nr=R.R.length;
  const rows=sp.map((s,i)=>{
    const id='f'+i,role=i<nr?'Reactant':'Product';
    let mid,hint='';
    if(S.mode==='worked'){mid=`<span class="big">${fx(s.f)}</span>`;hint=hintF(s);}
    else{
      mid=inp(id,s.f,{kind:'txt',wide:true,label:'Formula of '+s.name});
      if(S.mode==='guided')hint=hintF(s)+` <button class="btn sm" data-act="show" data-id="${id}">Show</button>`;
      else if(S.mode==='standard')hint=S.hintOpen[id]?hintF(s):`<button class="btn sm" data-act="hint" data-id="${id}">Hint</button>`;
    }
    return `<div class="frow"><div><b>${s.name}</b><div class="hint">${role}, ${({solid:'solid',liquid:'liquid',gas:'gas',aq:'in solution'})[s.form]}</div></div><div>${mid}</div><div class="hint">${hint}</div></div>`;
  }).join('');
  const ions=S.mode==='challenge'?'':`<details><summary>Ion charges you need to know</summary><div class="ionsets">
    <div class="ionset pos"><h3>Positive ions (cations)</h3><div class="ionrow"><b>1+</b><span>Na⁺ sodium</span><span>K⁺ potassium</span><span>Ag⁺ silver</span><span>NH₄⁺ ammonium</span><span>H⁺ hydrogen</span></div><div class="ionrow"><b>2+</b><span>Mg²⁺ magnesium</span><span>Ca²⁺ calcium</span><span>Zn²⁺ zinc</span><span>Cu²⁺ copper(II)</span><span>Fe²⁺ iron(II)</span></div><div class="ionrow"><b>3+</b><span>Al³⁺ aluminium</span><span>Fe³⁺ iron(III)</span></div></div>
    <div class="ionset neg"><h3>Negative ions (anions)</h3><div class="ionrow"><b>1−</b><span>Cl⁻ chloride</span><span>OH⁻ hydroxide</span><span>NO₃⁻ nitrate</span><span>CH₃COO⁻ ethanoate</span></div><div class="ionrow"><b>2−</b><span>O²⁻ oxide</span><span>S²⁻ sulfide</span><span>SO₄²⁻ sulfate</span><span>CO₃²⁻ carbonate</span></div><div class="ionrow"><b>3−</b><span>PO₄³⁻ phosphate</span></div></div></div></details>`;
  const btns=S.mode==='worked'?'':`<div class="chips">${S.mode!=='guided'?'<button class="btn teal" data-act="check" data-p="f">Check formulae</button>':''}<button class="btn" data-act="showall" data-p="f">Show all answers</button></div>`;
  $('stage').innerHTML=`<header><h2>Step 2 · Write the formulae</h2><p>${S.mode==='worked'?'Each formula below is built from the ions or the molecule. Read the reason beside it.':'Type the chemical formula for each substance. Subscripts can be typed as ordinary numbers, for example H2O or Ca(OH)2.'}</p></header>
  <div class="rows">${rows}</div>${btns}${ions}${navHTML()}`;
  paint();
}

/* ---------- step 3: balance ---------- */
function atomsHTML(){
  const vis=S.mode==='worked'||S.mode==='guided'||(S.mode==='standard'&&S.hintOpen.atoms)||S.checked.bal;
  const info=balInfo(),bad=info.find(x=>!x.ok);
  let msg='';
  if(vis||S.checked.bal){
    if(!bad)msg=simplest()?`<div class="callout good"><b>Balanced.</b> Every element has the same number of atoms on both sides, in the simplest ratio.</div>`:`<div class="callout warn"><b>Balanced, but not in the simplest ratio.</b> Divide every coefficient by the same number.</div>`;
    else if(S.mode==='guided'||S.mode==='standard')msg=`<div class="callout warn">${bad.e} is not balanced yet: ${bad.l} on the left and ${bad.r} on the right. Change a coefficient, never a subscript.</div>`;
    else msg=`<div class="callout bad">Not balanced yet.</div>`;
  }
  const tbl=vis?`<div class="grid"><table><thead><tr><th class="rl">Atoms of</th><th>Left</th><th>Right</th><th></th></tr></thead><tbody>${info.map(x=>`<tr><td class="rl"><b>${x.e}</b></td><td><span class="given">${x.l}</span></td><td><span class="given">${x.r}</span></td><td>${x.ok?'✓':'✗'}</td></tr>`).join('')}</tbody></table></div>`:'';
  return tbl+msg;
}
function stage3(){
  const R=S.R;
  if(!R.valid){$('stage').innerHTML=`<header><h2>Step 3 · Balance</h2></header><div class="callout bad">${R.msg}</div>${navHTML(false)}`;return;}
  const sp=allSp(),nr=R.R.length,w=S.mode==='worked';
  if(w)S.coefs=R.coefs.slice();
  const box=(s,i)=>`<div class="card" style="align-items:center;flex:0 0 auto"><div class="eq chem"><span>${fx(s.f)}</span></div>${w?`<span class="big">${S.coefs[i]}</span>`:`<div class="step"><button data-act="coef" data-i="${i}" data-d="-1" aria-label="Decrease coefficient of ${s.name}">−</button><output id="co${i}">${S.coefs[i]}</output><button data-act="coef" data-i="${i}" data-d="1" aria-label="Increase coefficient of ${s.name}">+</button></div>`}</div>`;
  const line=sp.map((s,i)=>(i===nr?'<span class="big" style="align-self:center">→</span>':i?'<span class="big plus" style="align-self:center">+</span>':'')+box(s,i)).join('');
  const method=w?`<div class="callout"><b>Method.</b> 1) Count the atoms of each element on both sides. 2) Pick an element that is out of balance, starting with one that appears in just one substance on each side. 3) Put a number in front of a formula to fix it. 4) Re-count everything, because that number changes other elements too. 5) Stop when every element matches in the smallest whole numbers.</div>`:'';
  const btns=w?'':`<div class="chips">${S.mode==='standard'?'<button class="btn" data-act="hint" data-id="atoms">Count the atoms</button>':''}${S.mode!=='guided'?'<button class="btn teal" data-act="checkbal">Check balance</button>':''}<button class="btn" data-act="solvebal">Show balanced equation</button></div>`;
  $('stage').innerHTML=`<header><h2>Step 3 · Balance the equation</h2><p>Use + and − to set the number in front of each formula. Atoms cannot be created or destroyed, so both sides must match.</p></header>
  <div class="chips" style="gap:10px;align-items:stretch">${line}</div>${method}${btns}<div id="atoms" class="rows">${atomsHTML()}</div>${navHTML()}`;
}
function refreshBal(){const a=$('atoms');if(a)a.innerHTML=atomsHTML();renderEq();renderRail();}

/* ---------- step 4: experiment bench ---------- */
const LIQ={CuSO4:'#3b9fe0','Cu(NO3)2':'#2f86c8',FeSO4:'#9fd49a'};
const liqColour=s=>{const f=s.f;return LIQ[f]||(s.ion&&s.ion.c.f==='Cu'?'#3b9fe0':s.ion&&s.ion.c.f==='Fe'&&s.ion.c.c===2?'#9fd49a':'#bfe3f0');};
const SOLID_COL={MgO:'white',ZnO:'white',CaO:'white',CuO:'black',Fe2O3:'red-brown',Al2O3:'white',Na2O:'white',K2O:'white',FeS:'black',CuS:'black',FeCl3:'dark red-brown',CuCl2:'yellow-brown',Cu:'pink-brown',Ag:'grey, shiny',Zn:'grey',Fe:'grey',Mg:'grey',C:'black'};
const solidCol=s=>SOLID_COL[s.f]||(s.ion&&s.ion.c.f==='Cu'?'blue-green':'white');
const SOLID_FILL={white:'#f3f3f3',black:'#2b2b2b','red-brown':'#9b3d22','pink-brown':'#c4673a','grey, shiny':'#b9bcc2',grey:'#9aa0a8','blue-green':'#4fb3b0','dark red-brown':'#6e2a1b','yellow-brown':'#b98a34'};
const solidFill=s=>SOLID_FILL[solidCol(s)]||'#f3f3f3';

function visHTML(i){
  const s=S.R.R[i],a=S.A[i];
  if(s.form==='aq'){
    const h=Math.min(1,a.V/250)*150,col=liqColour(s);
    let t='';for(let k=0;k<=5;k++)t+=`<line x1="30" x2="38" y1="${170-k*30}" y2="${170-k*30}" stroke="currentColor" opacity=".5"/><text x="26" y="${173-k*30}" font-size="9" text-anchor="end" fill="currentColor" opacity=".7">${k*50}</text>`;
    return `<svg viewBox="0 0 110 190" width="90" role="img" aria-label="Measuring cylinder holding ${a.V} cm³"><rect x="40" y="${170-h}" width="36" height="${h}" fill="${col}" opacity=".75"/><path d="M40 16V172H76V16" fill="none" stroke="currentColor" stroke-width="2"/>${t}<rect x="34" y="172" width="48" height="6" rx="2" fill="currentColor" opacity=".5"/></svg><div class="read">${fmt(a.V)} cm³</div><div class="hint">${fmt(a.c)} mol/dm³</div>`;
  }
  if(s.form==='gas'){
    const x=20+Math.min(1,a.V/20)*130;
    return `<svg viewBox="0 0 220 70" width="210" role="img" aria-label="Gas syringe holding ${a.V} dm³"><rect x="20" y="18" width="${x-20}" height="30" fill="var(--teal)" opacity=".25"/><rect x="20" y="18" width="140" height="30" fill="none" stroke="currentColor" stroke-width="2"/><rect x="${x-3}" y="12" width="6" height="42" fill="currentColor"/><line x1="${x}" x2="${x+50}" y1="33" y2="33" stroke="currentColor" stroke-width="4"/><rect x="${x+50}" y="22" width="6" height="22" fill="currentColor"/><line x1="160" x2="190" y1="33" y2="33" stroke="currentColor" stroke-width="2" opacity=".6"/></svg><div class="read">${fmt(a.V)} dm³</div><div class="hint">${fmt(a.V*1000)} cm³ at room temperature and pressure</div>`;
  }
  const hh=Math.min(1,a.m/40)*26+2,col=solidFill(s);
  return `<svg viewBox="0 0 180 90" width="170" role="img" aria-label="Balance reading ${a.m} g"><rect x="10" y="62" width="160" height="22" rx="6" fill="currentColor" opacity=".18"/><rect x="38" y="54" width="104" height="8" rx="3" fill="currentColor" opacity=".45"/><ellipse cx="90" cy="${54-hh/2}" rx="${14+hh*1.3}" ry="${hh/2+1}" fill="${col}" stroke="currentColor" stroke-opacity=".4"/><text x="150" y="78" font-size="11" text-anchor="end" font-family="monospace" fill="currentColor">${fmt(a.m)} g</text></svg><div class="read">${fmt(a.m)} g</div><div class="hint">balance reading${s.form==='liquid'?' (liquid poured into a pot)':''}</div>`;
}
function cardHTML(i){
  const s=S.R.R[i],a=S.A[i],aq=s.form==='aq',gas=s.form==='gas';
  const steps=gas?[[-1],[-0.1],[0.1],[1],[5]]:aq?[[-10],[-1],[1],[10]]:[[-1],[-0.1],[0.1],[1],[5]];
  const key=aq||gas?'V':'m';
  const unit=gas?'dm³':aq?'cm³':'g';
  const concs=[0.05,0.1,0.2,0.25,0.5,1,1.5,2,3];if(!concs.includes(a.c)&&aq)concs.push(a.c);concs.sort((x,y)=>x-y);
  const hold=aq?`<button class="btn sm teal" data-hold data-i="${i}" data-k="V" data-d="1" aria-label="Pour ${s.name}, hold to keep pouring">Pour (hold)</button>`:'';
  return `<div class="card"><h3><span>${fx(s.f)}</span> <span class="hint" style="font-weight:400">${s.name}</span></h3>
  <div id="vis${i}">${visHTML(i)}</div>
  ${aq?`<label class="f">Bottle label, concentration (mol/dm³)<select data-amt-c="${i}">${concs.map(c=>`<option value="${c}"${c===a.c?' selected':''}>${fmt(c)}</option>`).join('')}</select></label>`:''}
  <div class="btnrow">${hold}${steps.map(([d])=>`<button class="btn sm" data-act="adj" data-i="${i}" data-k="${key}" data-d="${d}">${d>0?'+':'−'}${Math.abs(d)}</button>`).join('')}<button class="btn sm" data-act="adj" data-i="${i}" data-k="${key}" data-d="-9999">Empty</button></div>
  <label class="f">Or type the amount (${unit})<input type="number" id="in${i}" data-amt="${i}" data-k="${key}" step="${gas?0.1:aq?1:0.1}" min="0" value="${a[key]}"></label></div>`;
}
function adj(i,k,d){
  const a=S.A[i],max=k==='V'?(S.R.R[i].form==='gas'?20:250):80;
  a[k]=Math.min(max,Math.max(0,Math.round((a[k]+d)*100)/100));
  $('vis'+i).innerHTML=visHTML(i);const el=$('in'+i);if(el)el.value=a[k];
  clearCalc();
}
function clearCalc(){
  Object.keys(S.vals).forEach(k=>{if(k[0]==='c')delete S.vals[k];});
  Object.keys(S.checked).forEach(k=>{if(k[0]==='c')delete S.checked[k];});
  Object.keys(S.shown).forEach(k=>{if(k[0]==='c')delete S.shown[k];});
  S.run=null;S.yld=null;S.reveal=0;const o=$('runout');if(o)o.innerHTML='';
}
const canTit=()=>S.cat==='neutralisation'&&S.R.valid&&S.R.R[0].form==='aq'&&S.R.R[1].form==='aq';

function stage4(){
  const R=S.R;
  if(!R.valid){$('stage').innerHTML=`<header><h2>Step 4 · Experiment</h2></header><div class="callout bad">${R.msg}</div>${navHTML(false)}`;return;}
  const tit=S.prac==='tit'&&canTit();
  const toggle=canTit()?`<div class="chips" role="group" aria-label="Practical"><button class="chip" data-act="prac" data-p="mix" aria-pressed="${!tit}">Mix and measure</button><button class="chip" data-act="prac" data-p="tit" aria-pressed="${tit}">Titration</button></div>`:'';
  const formulas=S.mode==='challenge'?'':`<div class="callout"><b>Amounts you can use.</b> <span class="tri"><span>solid or liquid</span><span>n = m ÷ Mr</span><span>solution</span><span>n = c × V ÷ 1000 (V in cm³)</span><span>gas</span><span>n = V ÷ 24 (V in dm³, room temperature and pressure)</span></span></div>`;
  let body;
  if(tit){body=titHTML();}
  else{
    body=`<div class="bench">${R.R.map((s,i)=>`<div id="card${i}">${cardHTML(i)}</div>`).join('')}</div>
    <div class="chips"><button class="btn" data-act="randamt">New random amounts</button><button class="btn pri" data-act="run">Mix and run the reaction</button></div>
    <div id="runout">${runHTML()}</div>`;
  }
  $('stage').innerHTML=`<header><h2>Step 4 · Set up the experiment</h2><p>${tit?'Add the acid from the burette to the alkali in the flask until the indicator changes colour. The reading at that point is the titre.':'Weigh, pour or measure each reactant. Change any amount you like. Then run the reaction and compare what happens with your equation.'}</p></header>
  ${toggle}${formulas}${body}${navHTML()}`;
}

/* running the reaction */
const GASOBS={H2:'Fizzing. The gas burns with a squeaky pop (hydrogen).',CO2:'Gas turns limewater milky (carbon dioxide).',O2:'Gas relights a glowing splint (oxygen).',NO2:'Brown fumes of nitrogen dioxide.',CO:'Colourless toxic gas that burns with a blue flame.',SO2:'Colourless gas with a choking smell that turns acidified dichromate green.'};
function obsList(R,c){
  const o=[],exo=R.cat==='combustion'||R.cat==='neutralisation'||R.cat==='displacement'||R.cat==='synthesis';
  const solids=R.P.filter(s=>s.form==='solid'),gases=R.P.filter(s=>s.form==='gas');
  if(R.cat==='synthesis'){o.push(`The reactants glow or burn. A ${solidCol(R.P[0])} solid forms: ${R.P[0].name}.`);if(R.R[0].f==='Mg'&&R.R[1].f==='O2')o.push('Very bright white flame.');}
  if(R.cat==='combustion'){o.push('Burns with a flame and gives out heat.');if(R.P.some(s=>s.f==='H2O'))o.push('Droplets of water condense on a cold surface.');if(R.R[0].f==='C'||R.R[0].f==='S')o.push('The solid glows as it burns.');}
  if(R.cat==='decomposition'){const r=R.R[0];o.push(`The solid is heated strongly${r.f==='CuCO3'?' and changes from green to black':r.f==='Cu(NO3)2'?' and changes from blue to black':''}.`);if(r.f==='H2O2')o.push('A little manganese(IV) oxide speeds it up. Bubbling is vigorous.');if(r.f==='H2O')o.push('Hydrogen collects at the negative electrode and oxygen at the positive. Hydrogen is twice the volume of oxygen.');}
  if(R.cat==='displacement'){
    const deposit=R.P.find(s=>s.form==='solid');
    if(deposit)o.push(`A ${solidCol(deposit)} solid (${deposit.name}) coats the metal.`);
    const cu=R.R[1].ion&&R.R[1].ion.c.f==='Cu',fe=R.R[1].ion&&R.R[1].ion.c.f==='Fe';
    if(cu&&c.left[1]<1e-9)o.push('The blue colour of the solution fades.');
    o.push('The temperature rises, so the reaction is exothermic.');
  }
  if(R.cat==='neutralisation'){if(R.type==='metal')o.push('The metal fizzes and gets smaller as it dissolves.');o.push('The temperature rises, so the reaction is exothermic.');if(R.R[1].form==='solid')o.push(c.left[1]>1e-9?`Some ${R.R[1].name} is left undissolved, which shows the acid has been used up.`:'The solid dissolves completely.');}
  gases.forEach(s=>{if(GASOBS[s.f])o.push(GASOBS[s.f]);});
  solids.filter(s=>R.cat!=='synthesis'&&R.cat!=='displacement').forEach(s=>o.push(`${solidCol(s)[0].toUpperCase()+solidCol(s).slice(1)} solid: ${s.name}.`));
  R.R.forEach((s,i)=>{if(c.left[i]>1e-6)o.push(`${s.name} is in excess, so some is left over.`);});
  if(c.allLim)o.push('Exactly the right amounts: nothing is left over.');
  return o;
}
function flaskSVG(R,c){
  const heat=['burn','heat'].includes(R.cond),aq=R.R.some(s=>s.form==='aq')||R.P.some(s=>s.form==='aq');
  const gas=R.P.some(s=>s.form==='gas'),sol=R.P.find(s=>s.form==='solid');
  let col='#bfe3f0';
  const cuLeft=R.R.some((s,i)=>s.ion&&s.ion.c.f==='Cu'&&c.left[i]>1e-9);
  const ap=R.P.find(s=>s.form==='aq');
  col=cuLeft?'#3b9fe0':ap?liqColour(ap):'#bfe3f0';
  return `<svg viewBox="0 0 170 220" width="170" role="img" aria-label="Reaction vessel after the reaction"><defs><clipPath id="fl"><path d="M64 14H106V74L146 154H24L64 74Z"/></clipPath></defs>
  ${aq?`<g clip-path="url(#fl)"><rect x="0" y="96" width="170" height="70" fill="${col}" opacity=".7"/></g>`:''}
  ${sol?`<g clip-path="url(#fl)"><rect x="0" y="144" width="170" height="12" fill="${solidFill(sol)}"/></g>`:''}
  <path d="M64 14H106V74L146 154H24L64 74Z" fill="none" stroke="currentColor" stroke-width="2.5"/>
  ${gas?[78,92,100,86].map((x,k)=>`<circle class="bubble" cx="${x+k*2}" cy="${140-k*4}" r="${3+k%2}" fill="none" stroke="currentColor" stroke-opacity=".7"/>`).join(''):''}
  ${heat?`<g class="flame"><path d="M85 214C70 200 78 184 85 170C92 184 100 200 85 214Z" fill="#f08a1c"/><path d="M85 212C78 204 82 194 85 188C88 194 92 204 85 212Z" fill="#ffd45a"/></g>`:''}</svg>`;
}
function runHTML(){
  const R=S.R;if(!R.valid||!S.run)return '';
  const c=S.run.c;
  if(c.bad)return `<div class="callout warn">Add some of every reactant first. At least one amount is zero.</div>`;
  const meas=[];
  R.P.forEach((s,j)=>{
    if(s.form==='gas')meas.push(`Gas syringe: <b>${fmt(c.nP[j]*24000)} cm³</b> of ${s.name}`);
    else if(s.form==='solid')meas.push(`Mass of ${s.name} collected: <b>${fmt(c.pAmt[j])} g</b>`);
  });
  R.R.forEach((s,i)=>{if(c.left[i]>1e-6&&s.form==='solid')meas.push(`Left over ${s.name}: <b>${fmt(c.leftAmt[i])} g</b>`);});
  if(R.cat==='synthesis'||R.cat==='decomposition'){const m0=S.A.reduce((t,a,i)=>t+(S.R.R[i].form==='solid'||S.R.R[i].form==='liquid'?a.m:0),0);meas.push(`Starting mass of solid on the balance: <b>${fmt(m0)} g</b>`);}
  return `<div class="vessel"><div>${flaskSVG(R,c)}</div><div style="flex:1;min-width:230px"><h3>What you see</h3><ul class="obs">${obsList(R,c).map(t=>`<li>${t}</li>`).join('')}</ul>${meas.length?`<h3 style="margin-top:12px">Measurements</h3><ul class="obs">${meas.map(t=>`<li>${t}</li>`).join('')}</ul>`:''}<p class="hint" style="margin-top:10px;color:var(--muted)">Now use these amounts in Step 5 to predict every number above before you check.</p></div></div>`;
}

/* ---------- titration ---------- */
function setupTit(){
  const R=S.R,co=R.coefs,keep=S.A[1].c,list=[0.05,0.1,0.2,0.25,0.5,1];
  const ve=(cb,ca)=>cb*25*co[0]/co[1]/ca;
  let best=null;
  const tryPair=(cb,ca)=>{const v=ve(cb,ca);if(v>=15&&v<=40&&(!best||Math.abs(v-25)<Math.abs(best.v-25)))best={cb,ca,v};};
  list.concat([2]).forEach(ca=>tryPair(keep,ca));
  if(!best)list.forEach(cb=>list.concat([2]).forEach(ca=>tryPair(cb,ca)));
  best=best||{cb:0.1,ca:0.1,v:ve(0.1,0.1)};
  S.A[1]={c:best.cb,V:25};S.A[0]={c:best.ca,V:0};
  const bt=BASES[R.R[1].f][1];
  S.tit={Ve:best.v,added:0,ind:(bt==='carb'||bt==='nh3')?'methyl':'phenol',done:false,titre:null,msg:''};
  clearCalc();
}
function indColour(){
  const t=S.tit,d=t.added-t.Ve;
  if(t.ind==='phenol'){if(d<-0.03)return `rgba(214,28,140,${(0.4+0.6*Math.max(0,1-t.added/t.Ve)).toFixed(2)})`;return Math.abs(d)<=0.03?'rgba(214,28,140,0.14)':'rgba(255,255,255,0.0)';}
  return d<-0.03?'#f2c200':Math.abs(d)<=0.03?'#f08a1c':'#e23b2a';
}
function titVis(){
  const t=S.tit,y=20+Math.min(1,t.added/50)*190;
  return `<svg viewBox="0 0 200 330" width="190" role="img" aria-label="Burette reading ${t.added.toFixed(2)} cm³ above a flask"><rect x="70" y="10" width="26" height="200" fill="none" stroke="currentColor" stroke-width="2"/><rect x="71" y="${y}" width="24" height="${210-y}" fill="#bfe3f0" opacity=".8"/><rect x="71" y="11" width="24" height="${Math.max(0,y-11)}" fill="none"/>
  ${[0,10,20,30,40,50].map(k=>`<line x1="96" x2="104" y1="${20+k*3.8}" y2="${20+k*3.8}" stroke="currentColor"/><text x="108" y="${23+k*3.8}" font-size="9" fill="currentColor">${k}</text>`).join('')}
  <path d="M83 210V222" stroke="currentColor" stroke-width="3"/>${t.added>0&&!t.done?'<circle cx="83" cy="228" r="2.5" fill="#bfe3f0"/>':''}
  <path d="M60 250H106V274L132 320H34L60 274Z" fill="${indColour()}" stroke="currentColor" stroke-width="2.5"/></svg>`;
}
function titHTML(){
  const R=S.R;if(!S.tit)setupTit();
  const t=S.tit,ac=R.R[0],bs=R.R[1];
  const hint=S.mode==='worked'?`The end-point is at about ${fmt(t.Ve)} cm³.`:S.mode==='guided'?`Hint: expect the colour to change somewhere between ${Math.floor(t.Ve-4)} and ${Math.ceil(t.Ve+4)} cm³. Go slowly near the end.`:'';
  return `<div class="cols"><div class="card" style="align-items:center"><div id="titvis">${titVis()}</div></div>
  <div style="flex:1;min-width:240px;display:flex;flex-direction:column;gap:10px">
   <div class="callout"><b>Burette:</b> ${ac.name}, ${fmt(S.A[0].c)} mol/dm³.<br><b>Flask:</b> 25.0 cm³ of ${bs.name}${S.mode==='challenge'||S.mode==='standard'?', concentration <b>unknown</b>':`, ${fmt(S.A[1].c)} mol/dm³`}, plus indicator.</div>
   <label class="f">Indicator<select data-tit-ind><option value="phenol"${t.ind==='phenol'?' selected':''}>Phenolphthalein (pink in alkali, colourless in acid)</option><option value="methyl"${t.ind==='methyl'?' selected':''}>Methyl orange (yellow in alkali, red in acid)</option></select></label>
   <div class="read" id="titread">Burette reading: ${t.added.toFixed(2)} cm³</div>
   <div class="btnrow chips"><button class="btn sm" data-act="titadd" data-d="5">+5</button><button class="btn sm" data-act="titadd" data-d="1">+1</button><button class="btn sm" data-act="titadd" data-d="0.1">+0.1</button><button class="btn sm" data-act="titadd" data-d="0.05">1 drop</button><button class="btn sm" data-act="titreset">Refill burette</button><button class="btn teal" data-act="titrecord">Record titre</button></div>
   ${hint?`<p class="hint" style="color:var(--muted)">${hint}</p>`:''}
   <div id="titmsg">${t.msg||''}</div></div></div>`;
}
function titRefresh(){
  const t=S.tit;$('titvis').innerHTML=titVis();$('titread').textContent=`Burette reading: ${t.added.toFixed(2)} cm³`;$('titmsg').innerHTML=t.msg||'';
}
function titAdd(d){
  const t=S.tit;if(t.done)return;
  t.added=Math.min(50,Math.round((t.added+d)*100)/100);
  t.msg=t.added>t.Ve+5?'<div class="callout warn">You have gone well past the end-point. Refill the burette and try again, adding drops near the end.</div>':'';
  titRefresh();
}
function titRecord(){
  const t=S.tit,d=t.added-t.Ve;
  if(t.added===0){t.msg='<div class="callout warn">Add some acid first.</div>';}
  else if(d<-0.1)t.msg='<div class="callout warn">The indicator has not reached its end-point colour yet. Keep adding acid.</div>';
  else if(d>0.1){t.msg='<div class="callout bad">Overshot: the indicator went past the end-point before you recorded. Refill the burette and repeat, using drops near the end.</div>';}
  else{t.done=true;t.titre=t.added;S.A[0].V=t.added;t.msg=`<div class="callout good"><b>Titre recorded: ${t.added.toFixed(2)} cm³.</b> Go to Step 5 to work out the unknown concentration.</div>`;renderRail();}
  titRefresh();
}

/* ---------- step 5: the equation grid ---------- */
function mrWork(f){
  const p=parse(f);return Object.keys(p).map(e=>`${p[e]>1?p[e]+'×':''}${AR[e]}`).join(' + ')+' = '+fmt(mr(f));
}
function rowsDef(c){
  const R=S.R,sp=allSp(),nr=R.R.length,co=R.coefs,A=S.A,isR=i=>i<nr,jj=i=>i-nr;
  const FN={m:'n = m ÷ Mr',cv:'n = c × V ÷ 1000',v:'n = V ÷ 24'};
  const fn=i=>({solid:'m',liquid:'m',gas:'v',aq:'cv'})[sp[i].form];
  const nTxt=i=>{const s=sp[i],a=A[i];return s.form==='aq'?`${fmt(a.c)} × ${fmt(a.V)} ÷ 1000`:s.form==='gas'?`${fmt(a.V)} ÷ 24`:`${fmt(a.m)} ÷ ${fmt(mr(s.f))}`;};
  const unit=i=>sp[i].form==='gas'?'dm³':sp[i].form==='aq'?'g of solute':'g';
  return [
   {k:'co',label:'Coefficient',why:'The balanced equation gives the mole ratio.',hint:'Read the number in front of each formula in your balanced equation.',app:()=>true,ans:i=>co[i],tol:.001,blank:m=>m==='challenge',text:i=>co[i]},
   {k:'given',label:'Amount given',why:'These are the amounts you set in Step 4.',app:isR,text:i=>{const s=sp[i],a=A[i];return s.form==='aq'?`${fmt(a.c)} mol/dm³ × ${fmt(a.V)} cm³`:s.form==='gas'?`${fmt(a.V)} dm³`:`${fmt(a.m)} g`;},fixed:true},
   {k:'mr',label:'Mr',why:'Add the Ar of every atom in the formula.',hint:'Add up the Ar of every atom. Brackets multiply everything inside them.',app:()=>true,ans:i=>mr(sp[i].f),blank:m=>m!=='guided',text:i=>fmt(mr(sp[i].f)),work:i=>mrWork(sp[i].f)},
   {k:'fn',label:'Which equation?',sel:true,app:isR,ans:i=>fn(i),opts:[['','choose…'],['m',FN.m],['cv',FN.cv],['v',FN.v]],blank:()=>true},
   {k:'n',label:'Moles available',why:'Change each reactant amount into moles.',hint:'Solids and liquids: n = m ÷ Mr. Solutions: n = c × V ÷ 1000 (V in cm³). Gases: n = V ÷ 24 (V in dm³).',app:isR,ans:i=>c.n[i],blank:()=>true,text:i=>fmt(c.n[i]),work:i=>nTxt(i)+' = '+fmt(c.n[i]),cap:i=>FN[fn(i)]},
   {k:'ratio',label:'Moles ÷ coefficient',why:'Dividing by the coefficient makes the reactants comparable. The smallest value runs out first.',hint:'Divide each reactant’s moles by its own coefficient.',app:isR,ans:i=>c.ratio[i],blank:()=>true,text:i=>fmt(c.ratio[i]),work:i=>`${fmt(c.n[i])} ÷ ${co[i]} = ${fmt(c.ratio[i])}`,cap:()=>'moles ÷ coefficient'},
   {k:'lim',label:'Limiting or in excess?',sel:true,why:'The reactant with the smallest ratio is limiting. It controls how much product forms.',hint:'Smallest (moles ÷ coefficient) means limiting. If they are equal, both are limiting.',app:isR,ans:i=>c.lim[i]?'lim':'exc',opts:[['','choose…'],['lim','Limiting'],['exc','In excess']],blank:()=>true,text:i=>c.lim[i]?'Limiting':'In excess'},
   {k:'used',label:'Moles that react or form',why:'Multiply each coefficient by the smallest ratio. That is how many moles of each substance take part.',hint:'moles = coefficient × the smallest (moles ÷ coefficient).',app:()=>true,ans:i=>isR(i)?c.used[i]:c.nP[jj(i)],blank:()=>true,text:i=>fmt(isR(i)?c.used[i]:c.nP[jj(i)]),work:i=>`${co[i]} × ${fmt(c.xi)} = ${fmt(isR(i)?c.used[i]:c.nP[jj(i)])}`,cap:()=>'coefficient × smallest ratio'},
   {k:'left',label:'Moles left over',why:'Moles available minus moles used. The limiting reactant has none left.',hint:'moles left = moles available − moles that react.',app:isR,ans:i=>c.left[i],blank:()=>true,text:i=>fmt(c.left[i]),work:i=>`${fmt(c.n[i])} − ${fmt(c.used[i])} = ${fmt(c.left[i])}`,cap:()=>'available − used'},
   {k:'amt',label:'Left over (reactants) or made (products)',why:'Convert moles back to a mass or a gas volume.',hint:'mass = n × Mr, gas volume = n × 24 dm³. For a dissolved substance this is the mass of solute.',app:()=>true,ans:i=>isR(i)?c.leftAmt[i]:c.pAmt[jj(i)],unit,blank:()=>true,text:i=>fmt(isR(i)?c.leftAmt[i]:c.pAmt[jj(i)])+' '+unit(i),work:i=>{const n=isR(i)?c.left[i]:c.nP[jj(i)];return sp[i].form==='gas'?`${fmt(n)} × 24 = ${fmt(n*24)} dm³`:`${fmt(n)} × ${fmt(mr(sp[i].f))} = ${fmt(n*mr(sp[i].f))} g`;},cap:i=>sp[i].form==='gas'?'n × 24':'n × Mr'},
   {k:'conc',label:'Concentration of the solution made',why:'Concentration is moles ÷ volume. Assume the solution volumes simply add.',hint:'c = n ÷ (total cm³ ÷ 1000). The total volume is the sum of the solution volumes you used.',app:i=>!isR(i)&&sp[i].form==='aq'&&c.pConc[jj(i)]!=null,ans:i=>c.pConc[jj(i)],blank:()=>true,text:i=>fmt(c.pConc[jj(i)])+' mol/dm³',work:i=>`${fmt(c.nP[jj(i)])} ÷ (${fmt(c.Vt)} ÷ 1000) = ${fmt(c.pConc[jj(i)])}`,unit:()=>'mol/dm³',cap:()=>'n ÷ V in dm³'}
  ];
}
function stage5(){
  const R=S.R;
  if(!R.valid){$('stage').innerHTML=`<header><h2>Step 5 · Calculate</h2></header><div class="callout bad">${R.msg}</div>${navHTML(false)}`;return;}
  if(S.prac==='tit'&&canTit())return stage5tit();
  const c=calc(R,S.A),m=S.mode,sp=allSp(),nr=R.R.length;
  if(c.bad){$('stage').innerHTML=`<header><h2>Step 5 · Calculate</h2></header><div class="callout warn">Every reactant needs an amount above zero. Go back to Step 4 and add some.</div>${navHTML()}`;return;}
  let rows=rowsDef(c).filter(r=>(r.k!=='fn'||m==='challenge')&&sp.some((s,i)=>r.app(i)));
  S.rowMeta=[];
  const total=rows.length,shownRows=m==='worked'?Math.min(total,S.reveal+1):total;
  const head=`<tr><th class="rl">Working</th>${sp.map((s,i)=>(i===nr?'<th class="arrowcol">→</th>':'')+`<th>${fx(s.f)}</th>`).join('')}</tr>`;
  const body=rows.map((r,ri)=>{
    const ids=[];
    const cells=sp.map((s,i)=>{
      const pre=i===nr?'<td class="arrowcol"></td>':'';
      if(!r.app(i))return pre+'<td>–</td>';
      const id=`c_${r.k}_${i}`,w=m==='worked';
      if(r.fixed||w||(r.blank&&!r.blank(m))){
        if(r.sel&&w)return pre+`<td><span class="tag ${c.lim[i]?'lim':'exc'}">${r.text(i)}</span></td>`;
        return pre+`<td><span class="given">${r.text(i)}</span>${w&&r.work?`<small style="display:block;color:var(--muted);font-size:11.5px">${r.work(i)}</small>`:''}</td>`;
      }
      ids.push(id);
      if(r.sel)return pre+`<td>${selF(id,r.ans(i),r.opts,{label:r.label+' for '+s.name})}</td>`;
      return pre+`<td>${inp(id,r.ans(i),{tol:r.tol,label:r.label+' for '+s.name})}${r.unit?`<small style="display:block;color:var(--muted)">${r.unit(i)}</small>`:''}${m==='guided'&&r.cap?`<small style="display:block;color:var(--muted)">${r.cap(i)}</small>`:''}</td>`;
    }).join('');
    S.rowMeta.push({k:r.k,ids});
    let small='';
    if(m==='worked')small=r.why||'';
    else if(m==='guided')small=r.hint||r.why||'';
    else if(m==='standard'&&r.hint&&ids.length)small=S.hintOpen['r'+r.k]?r.hint:`<button class="btn sm" data-act="hint" data-id="r${r.k}">Hint</button>`;
    let btn='';
    if(ids.length&&m==='standard')btn=`<div class="chips" style="margin-top:4px"><button class="btn sm" data-act="checkrow" data-r="${r.k}">Check</button><button class="btn sm" data-act="showrow" data-r="${r.k}">Show</button></div>`;
    if(ids.length&&m==='guided')btn=`<div style="margin-top:4px"><button class="btn sm" data-act="showrow" data-r="${r.k}">Show row</button></div>`;
    return `<tr data-row="${r.k}"${ri>=shownRows?' hidden':''}><td class="rl"><b>${r.label}</b><small>${small}</small>${btn}</td>${cells}</tr>`;
  }).join('');
  const limNames=R.R.filter((s,i)=>c.lim[i]).map(s=>fx(s.f)).join(' and ');
  const exNames=R.R.filter((s,i)=>!c.lim[i]).map(s=>fx(s.f)).join(' and ');
  const summary=`<div class="callout good"><b>Result.</b> ${c.allLim?'The reactants are in exactly the right ratio, so none is in excess.':`${limNames} is the limiting reagent and ${exNames} is in excess.`} The reaction makes ${R.P.map((s,j)=>`${fmt(c.pAmt[j])} ${unitOf(s)==='dm³'?'dm³':'g'} of ${fx(s.f)}`).join(', ')}.</div>`;
  const wbtn=m==='worked'?`<div class="chips">${S.reveal+1<total?`<button class="btn pri" data-act="revnext">Show next step ▸</button><button class="btn" data-act="revall">Show every step</button>`:''}<button class="btn" data-act="revreset">Start again</button></div>`:'';
  const gbtn=m==='worked'?'':`<div class="chips">${m!=='guided'?'<button class="btn teal" data-act="check" data-p="c_">Check all</button>':''}<button class="btn" data-act="showall" data-p="c_">Show all answers</button><button class="btn" data-act="randamt">New amounts</button></div>`;
  const intro={worked:'Read down the columns. Each step has its working and the reason beside it.',guided:'Each row unlocks when the row above is right. Answers within 2% count, so round sensibly.',standard:'Fill in the blanks. Check a row when you are ready, or check everything at the end.',challenge:'Choose the equation, then work out every blank. Press Check all when you are finished.'}[m];
  $('stage').innerHTML=`<header><h2>Step 5 · Work it out under the equation</h2><p>${intro}</p></header>
  <div class="grid"><table><thead>${head}</thead><tbody>${body}</tbody></table></div>
  ${wbtn}${gbtn}${(m==='worked'&&S.reveal+1>=total)?summary:'<div id="final"></div>'}${yieldHTML(c)}
  ${S.run&&!S.run.c.bad?`<details open><summary>Bench results from Step 4</summary><div style="margin-top:10px">${runHTML()}</div></details>`:''}${navHTML()}`;
  paint();
}
function updateLocks(){
  if(S.step!==5||S.mode!=='guided'||!S.rowMeta)return;
  let prevDone=true;
  S.rowMeta.forEach(rm=>{
    const tr=document.querySelector(`tr[data-row="${rm.k}"]`);if(!tr)return;
    const locked=!prevDone&&rm.ids.length>0;
    tr.classList.toggle('locked',locked);
    tr.querySelectorAll('[data-f]').forEach(el=>{el.disabled=locked;});
    if(rm.ids.length)prevDone=prevDone&&rm.ids.every(id=>isRight(id)===true||S.shown[id]);
  });
}
function stage5tit(){
  const R=S.R,m=S.mode,co=R.coefs;
  if(!S.tit||!S.tit.done){$('stage').innerHTML=`<header><h2>Step 5 · Titration calculation</h2></header><div class="callout warn">Finish the titration and record the titre in Step 4 first.</div>${navHTML()}`;return;}
  const ca=S.A[0].c,Va=S.A[0].V,Vb=25,na=ca*Va/1000,nb=na*co[1]/co[0],cb=nb/(Vb/1000);
  const steps=[
   {k:'t1',label:'Moles of acid used',why:'n = c × V ÷ 1000, using the titre as V in cm³.',ans:na,work:`${fmt(ca)} × ${fmt(Va)} ÷ 1000 = ${fmt(na)}`,unit:'mol'},
   {k:'t2',label:'Moles of alkali that reacted',why:`The equation says ${co[0]} mol of acid reacts with ${co[1]} mol of ${R.R[1].name}. Scale the moles of acid by this ratio.`,ans:nb,work:`${fmt(na)} × ${co[1]} ÷ ${co[0]} = ${fmt(nb)}`,unit:'mol'},
   {k:'t3',label:'Concentration of the alkali',why:'c = n ÷ (V ÷ 1000), with V = 25.0 cm³.',ans:cb,work:`${fmt(nb)} ÷ (25.0 ÷ 1000) = ${fmt(cb)}`,unit:'mol/dm³'}
  ];
  const rows=steps.map(st=>{
    let mid;
    if(m==='worked')mid=`<span class="given">${st.work}</span>`;
    else mid=inp('c_'+st.k,st.ans,{label:st.label})+` <span class="hint">${st.unit}</span>`;
    const hint=m==='worked'||m==='guided'?st.why:m==='standard'?(S.hintOpen['r'+st.k]?st.why:`<button class="btn sm" data-act="hint" data-id="r${st.k}">Hint</button>`):'';
    return `<div class="frow"><div><b>${st.label}</b></div><div>${mid}</div><div class="hint">${hint}</div></div>`;
  }).join('');
  const btns=m==='worked'?'':`<div class="chips">${m!=='guided'?'<button class="btn teal" data-act="check" data-p="c_t">Check</button>':''}<button class="btn" data-act="showall" data-p="c_t">Show answers</button></div>`;
  $('stage').innerHTML=`<header><h2>Step 5 · Work out the unknown concentration</h2><p>Your titre was ${fmt(Va)} cm³ of ${fmt(ca)} mol/dm³ ${R.R[0].name} reacting with 25.0 cm³ of ${R.R[1].name}.</p></header>
  <div class="callout">Equation: ${R.R.map((s,i)=>`${R.coefs[i]>1?R.coefs[i]+' ':''}${fx(s.f)}`).join(' + ')} → ${R.P.map((s,j)=>`${R.coefs[R.R.length+j]>1?R.coefs[R.R.length+j]+' ':''}${fx(s.f)}`).join(' + ')}</div>
  <div class="rows">${rows}</div>${btns}<div class="callout good">The bottle was labelled ${fmt(S.A[1].c)} mol/dm³. Titres are only accurate to about ±0.1 cm³, so a small difference is normal.</div>${navHTML()}`;
  paint();
}

/* ---------- mole practice ---------- */
const T={kind:'mass',q:null,val:'',res:'',hint:false,shown:false,ok:0,n:0,done:false};
const POOL=['NaCl','CaCO3','H2SO4','HCl','NaOH','MgO','CuSO4','Al2O3','H2O','CO2','KNO3','NH3','Mg(OH)2','Na2CO3','Ca(OH)2','HNO3'];
const TRI={mass:'n = m ÷ Mr   ·   m = n × Mr   ·   Mr = m ÷ n',conc:'n = c × V ÷ 1000 (V in cm³)   ·   c = n ÷ (V ÷ 1000)   ·   V = n ÷ c × 1000 (cm³)',gas:'n = V ÷ 24 (V in dm³)   ·   V = n × 24 (dm³)   ·   1 dm³ = 1000 cm³'};
function newQ(){
  const k=T.kind,v=Math.floor(Math.random()*(k==='gas'?4:3));let q;
  if(k==='mass'){
    const f=rnd(POOL),M=mr(f),n=rnd([0.05,0.1,0.2,0.25,0.5,1,1.5,2,3]),m=Math.round(n*M*100)/100,nn=m/M;
    q=[{t:`How many moles are in ${fmt(m)} g of ${fx(f)}?`,a:nn,u:'mol',w:`Mr = ${mrWork(f)}. n = m ÷ Mr = ${fmt(m)} ÷ ${fmt(M)} = ${fmt(nn)} mol`},
       {t:`What is the mass of ${fmt(n)} mol of ${fx(f)}?`,a:n*M,u:'g',w:`Mr = ${mrWork(f)}. m = n × Mr = ${fmt(n)} × ${fmt(M)} = ${fmt(n*M)} g`},
       {t:`${fmt(m)} g of a compound is ${fmt(n)} mol. What is its Mr?`,a:m/n,u:'',w:`Mr = m ÷ n = ${fmt(m)} ÷ ${fmt(n)} = ${fmt(m/n)}`}][v];
  }else if(k==='conc'){
    const c=rnd([0.1,0.2,0.5,1,2]),V=rnd([10,20,25,50,100,250]),n=c*V/1000;
    q=[{t:`How many moles are in ${V} cm³ of a ${fmt(c)} mol/dm³ solution?`,a:n,u:'mol',w:`n = c × V ÷ 1000 = ${fmt(c)} × ${V} ÷ 1000 = ${fmt(n)} mol`},
       {t:`${fmt(n)} mol of solute is dissolved to make ${V} cm³ of solution. What is the concentration?`,a:c,u:'mol/dm³',w:`c = n ÷ (V ÷ 1000) = ${fmt(n)} ÷ ${fmt(V/1000)} = ${fmt(c)} mol/dm³`},
       {t:`What volume of ${fmt(c)} mol/dm³ solution contains ${fmt(n)} mol? Answer in cm³.`,a:V,u:'cm³',w:`V = n ÷ c × 1000 = ${fmt(n)} ÷ ${fmt(c)} × 1000 = ${fmt(V)} cm³`}][v];
  }else{
    const n=rnd([0.05,0.1,0.25,0.5,1,2,3]),V=n*24;
    q=[{t:`How many moles are in ${fmt(V)} dm³ of gas at room temperature and pressure?`,a:n,u:'mol',w:`n = V ÷ 24 = ${fmt(V)} ÷ 24 = ${fmt(n)} mol`},
       {t:`What volume does ${fmt(n)} mol of gas occupy at room temperature and pressure? Answer in dm³.`,a:V,u:'dm³',w:`V = n × 24 = ${fmt(n)} × 24 = ${fmt(V)} dm³`},
       {t:`How many moles are in ${fmt(V*1000)} cm³ of gas at room temperature and pressure?`,a:n,u:'mol',w:`${fmt(V*1000)} cm³ = ${fmt(V)} dm³. n = ${fmt(V)} ÷ 24 = ${fmt(n)} mol`},
       {t:`What volume does ${fmt(n)} mol of gas occupy at room temperature and pressure? Answer in cm³.`,a:V*1000,u:'cm³',w:`V = n × 24 = ${fmt(V)} dm³ = ${fmt(V*1000)} cm³`}][v];
  }
  T.q=q;T.val='';T.res='';T.hint=false;T.shown=false;T.done=false;
}
function renderTool(){
  if(!T.q)newQ();
  const m=S.mode,q=T.q;
  const kinds=[['mass','Mass and Mr'],['conc','Concentration and volume'],['gas','Gas volume']];
  const showTri=m==='worked'||m==='guided'||T.hint;
  $('toolView').innerHTML=`<section class="stage"><header><h2>Mole practice</h2><p>Three routes to the same unit. Pick one and work through random questions. Answers within 2% are accepted.</p></header>
  <div class="chips" role="group" aria-label="Calculation type">${kinds.map(([k,l])=>`<button class="chip" data-act="tkind" data-k="${k}" aria-pressed="${T.kind===k}">${l}</button>`).join('')}</div>
  <div class="q"><div class="text">${q.t}</div>
   ${showTri?`<div class="callout"><span class="tri"><span>${TRI[T.kind]}</span></span></div>`:''}
   <div class="chips" style="align-items:center"><input class="blank" id="tin" data-tin type="text" inputmode="decimal" autocomplete="off" value="${esc(T.val)}" aria-label="Your answer"> <span>${q.u}</span>
    <button class="btn teal" data-act="tcheck">Check</button>
    ${m==='standard'&&!T.hint?'<button class="btn" data-act="thint">Hint</button>':''}
    <button class="btn" data-act="tshow">Show working</button><button class="btn pri" data-act="tnew">New question</button></div>
   <div id="tres">${toolRes()}</div></div>
  <div class="score">${T.ok} correct out of ${T.n} attempts</div></section>`;
}
function toolRes(){
  const q=T.q;let h=T.res||'';
  if(T.shown||S.mode==='worked')h+=`<div class="callout"><b>Working.</b> ${q.w}</div>`;
  return h;
}
function toolCheck(){
  const x=parseFloat(String(T.val).replace(',','.')),q=T.q;
  if(isNaN(x)){T.res='<div class="callout warn">Type a number first.</div>';return;}
  T.n++;
  if(Math.abs(x-q.a)<=Math.abs(q.a)*0.02+1e-9){if(!T.done)T.ok++;T.done=true;T.res='<div class="callout good"><b>Correct.</b></div>';}
  else T.res=`<div class="callout bad">Not quite. ${S.mode==='challenge'?'Try again.':'Check your formula and your units, then try again.'}</div>`;
}

/* ---------- balancing practice ---------- */
const B={filter:'all',level:'easy',R:null,co:[],checked:false,shown:false,hint:false,ok:0,n:0,done:false};
const BF=[['all','All types'],['synthesis','Synthesis'],['combustion','Combustion'],['decomposition','Decomposition'],['displacement','Displacement'],['neutralisation','Acid reactions']];
const BL=[['easy','Easy'],['medium','Medium'],['hard','Hard']];
function randPick(cat){
  const p={};
  if(cat==='neutralisation'){p.type=rnd(Object.keys(NTYPE));p.acid=rnd(OPTS.neutralisation[0].list)[0];p.base=rnd(NBASE[p.type]);}
  else OPTS[cat].forEach(o=>{p[o.k]=rnd(o.list)[0];});
  return p;
}
const eqKey=R=>R.R.concat(R.P).map(s=>s.f).join();
function newBal(){
  const fits=(R,relax)=>{const mx=Math.max(...R.coefs),ns=R.R.length+R.P.length;return relax||(B.level==='easy'?mx<=2:B.level==='medium'?(mx>=2&&mx<=4):(mx>=4||ns>=5));};
  for(let pass=0;pass<2;pass++)for(let k=0;k<250;k++){
    const cat=B.filter==='all'?rnd(Object.keys(CATS)):B.filter;
    const R=build(cat,randPick(cat));
    if(!R.valid||!fits(R,pass===1)||(B.R&&eqKey(B.R)===eqKey(R)))continue;
    B.R=R;B.co=R.R.concat(R.P).map(()=>1);B.checked=false;B.shown=false;B.hint=false;B.done=false;return;
  }
}
function balInfoB(){
  const R=B.R,sp=R.R.concat(R.P),nr=R.R.length,els=[...new Set(sp.flatMap(s=>Object.keys(parse(s.f))))];
  return els.map(e=>{let l=0,r=0;sp.forEach((s,i)=>{const c=(parse(s.f)[e]||0)*B.co[i];if(i<nr)l+=c;else r+=c;});return {e,l,r,ok:l===r};});
}
function renderBal(){
  if(!B.R)newBal();
  const R=B.R,sp=R.R.concat(R.P),nr=R.R.length,m=S.mode,w=m==='worked';
  if(w)B.co=R.coefs.slice();
  const info=balInfoB(),bad=info.find(x=>!x.ok),bal=!bad,simp=B.co.reduce(gcd)===1;
  const vis=w||m==='guided'||(m==='standard'&&B.hint)||B.checked||B.shown;
  let fb='';
  if(bal&&simp&&(w||B.checked||B.shown||m==='guided'))fb=`<div class="callout good"><b>Balanced.</b> Every element matches on both sides in the simplest ratio.</div>`;
  else if(bal&&(w||B.checked||B.shown||m==='guided'))fb=`<div class="callout warn"><b>Balanced, but not simplest.</b> Divide every coefficient by the same number.</div>`;
  else if(!bal&&(B.checked||(m==='guided')))fb=m==='challenge'?`<div class="callout bad">Not balanced yet.</div>`:`<div class="callout warn">${bad.e} is not balanced: ${bad.l} on the left, ${bad.r} on the right. Change a coefficient, never a subscript.</div>`;
  const box=(s,i)=>`<div class="card" style="align-items:center;flex:0 0 auto"><div class="eq chem"><span>${fx(s.f)}</span></div>${w?`<span class="big">${B.co[i]}</span>`:`<div class="step"><button data-act="bco" data-i="${i}" data-d="-1" aria-label="Decrease coefficient of ${s.name}">−</button><output>${B.co[i]}</output><button data-act="bco" data-i="${i}" data-d="1" aria-label="Increase coefficient of ${s.name}">+</button></div>`}</div>`;
  const line=sp.map((s,i)=>(i===nr?'<span class="big" style="align-self:center">→</span>':i?'<span class="big plus" style="align-self:center">+</span>':'')+box(s,i)).join('');
  const tbl=vis?`<div class="grid"><table><thead><tr><th class="rl">Atoms of</th><th>Left</th><th>Right</th><th></th></tr></thead><tbody>${info.map(x=>`<tr><td class="rl"><b>${x.e}</b></td><td><span class="given">${x.l}</span></td><td><span class="given">${x.r}</span></td><td>${x.ok?'✓':'✗'}</td></tr>`).join('')}</tbody></table></div>`:'';
  const method=w?`<div class="callout"><b>Method.</b> Count the atoms of each element on both sides. Fix one element at a time by changing a coefficient, starting with an element that appears in only one substance on each side. Re-count after every change. Finish by checking the numbers are the smallest whole numbers.</div>`:'';
  const btns=w?'':`<div class="chips"><button class="btn teal" data-act="bcheck">Check</button>${m==='standard'&&!B.hint?'<button class="btn" data-act="bhint">Count the atoms</button>':''}<button class="btn" data-act="bshow">Show answer</button></div>`;
  $('balView').innerHTML=`<section class="stage"><header><h2>Balancing practice</h2><p>Set the coefficient in front of each formula so both sides have the same atoms. Use the type and level buttons to choose what you practise.</p></header>
  <div class="chips" role="group" aria-label="Reaction type">${BF.map(([k,l])=>`<button class="chip" data-act="bfilter" data-k="${k}" aria-pressed="${B.filter===k}">${l}</button>`).join('')}</div>
  <div class="chips" role="group" aria-label="Level">${BL.map(([k,l])=>`<button class="chip" data-act="blevel" data-k="${k}" aria-pressed="${B.level===k}">${l}</button>`).join('')}</div>
  <div class="callout"><b>${CATS[R.cat].label}:</b> ${R.R.map(s=>s.name).join(' + ')} → ${R.P.map(s=>s.name).join(' + ')}</div>
  <div class="chips" style="gap:10px;align-items:stretch">${line}</div>${method}${btns}${fb}${tbl}
  <div class="nav"><span class="score">${B.ok} balanced correctly out of ${B.n} checks</span><span><button class="btn pri" data-act="bnew">New equation</button></span></div></section>`;
}

/* ---------- percentage yield in step 5 ---------- */
function yieldHTML(c){
  const R=S.R,key=JSON.stringify(S.A)+R.coefs.join(),m=S.mode;
  if(!S.yld||S.yld.key!==key){
    let j=R.P.findIndex(s=>s.form==='solid');if(j<0)j=R.P.findIndex(s=>s.form==='gas');if(j<0)j=R.P.findIndex(s=>s.form==='liquid');if(j<0)j=0;
    const th=c.pAmt[j],pct=rnd([48,55,62,68,74,80,86,92]),actual=parseFloat((th*pct/100).toPrecision(3));
    S.yld={key,j,th,actual,pct:actual/th*100};
  }
  const y=S.yld,s=R.P[y.j],u=s.form==='gas'?'dm³':'g';
  let mid,hint='';
  if(m==='worked'){mid=`<span class="given">${fmt(y.actual)} ÷ ${fmt(y.th)} × 100 = ${fmt(y.pct)} %</span>`;hint='Percentage yield = actual yield ÷ theoretical yield × 100.';}
  else{
    mid=inp('c_y1',y.pct,{label:'Percentage yield'})+' <span class="hint">%</span>';
    if(m==='guided')hint='Percentage yield = actual yield ÷ theoretical yield × 100. Both amounts must be in the same unit.';
    else if(m==='standard')hint=S.hintOpen.ry?'Percentage yield = actual yield ÷ theoretical yield × 100.':'<button class="btn sm" data-act="hint" data-id="ry">Hint</button>';
  }
  const btn=m==='worked'?'':`<div class="chips">${m!=='guided'?'<button class="btn teal" data-act="check" data-p="c_y">Check</button>':''}<button class="btn" data-act="showall" data-p="c_y">Show answer</button></div>`;
  return `<div class="card" style="margin-top:4px"><h3>Percentage yield</h3>
  <p>The ${fx(s.f)} you worked out above is the <b>theoretical yield</b>: ${fmt(y.th)} ${u}. In the experiment only <b>${fmt(y.actual)} ${u}</b> was collected. This is the <b>actual yield</b>. It is lower because some product is lost in transfers, in side reactions, or because the reaction does not go to completion.</p>
  <div class="frow"><div><b>Percentage yield</b></div><div>${mid}</div><div class="hint">${hint}</div></div>${btn}</div>`;
}

/* ---------- percentage yield practice ---------- */
const Y={kind:'simple',q:null,vals:{},checked:false,shown:false,hint:false,ok:0,n:0,done:false};
function newY(){
  const pcts=[40,45,50,55,60,65,70,75,80,85,90,95];
  if(Y.kind==='simple'){
    const theo=rnd([2.5,4,5,8,10,12.5,20,25,40,50]),pct=rnd(pcts),act=Math.round(theo*pct)/100,v=Math.floor(Math.random()*3);
    const pp=act/theo*100;
    Y.q=[{t:`A reaction should give ${fmt(theo)} g of product (the theoretical yield). Only ${fmt(act)} g is collected. What is the percentage yield?`,b:[{l:'Percentage yield',a:pp,u:'%'}],w:`% yield = actual ÷ theoretical × 100 = ${fmt(act)} ÷ ${fmt(theo)} × 100 = ${fmt(pp)} %`},
         {t:`The theoretical yield of a reaction is ${fmt(theo)} g. The percentage yield is ${pct}%. What mass of product is actually collected?`,b:[{l:'Actual yield',a:theo*pct/100,u:'g'}],w:`actual = theoretical × % yield ÷ 100 = ${fmt(theo)} × ${pct} ÷ 100 = ${fmt(theo*pct/100)} g`},
         {t:`A student collects ${fmt(act)} g of product. The percentage yield is ${pct}%. What is the theoretical yield?`,b:[{l:'Theoretical yield',a:act/pct*100,u:'g'}],w:`theoretical = actual ÷ % yield × 100 = ${fmt(act)} ÷ ${pct} × 100 = ${fmt(act/pct*100)} g`}][v];
  }else{
    for(let k=0;k<300;k++){
      const cat=rnd(Object.keys(CATS)),R=build(cat,randPick(cat));
      if(!R.valid||!(R.R[0].form==='solid'||R.R[0].form==='liquid'))continue;
      const js=R.P.map((s,j)=>j).filter(j=>R.P[j].form!=='aq');if(!js.length)continue;
      const j=rnd(js),r0=R.R[0],P=R.P[j],M=mr(r0.f),n=rnd([0.05,0.1,0.2,0.25,0.5]),mass=Math.round(n*M*10)/10;
      if(mass<0.5||mass>60)continue;
      const nr=R.R.length,nR=mass/M,nP=nR*R.coefs[nr+j]/R.coefs[0],theo=amtOf(P,nP),pct=rnd(pcts),act=parseFloat((theo*pct/100).toPrecision(3)),pp=act/theo*100;
      const u=P.form==='gas'?'dm³':'g';
      const eq=R.R.map((s,i)=>(R.coefs[i]>1?R.coefs[i]+' ':'')+fx(s.f)).join(' + ')+' → '+R.P.map((s,jj)=>(R.coefs[nr+jj]>1?R.coefs[nr+jj]+' ':'')+fx(s.f)).join(' + ');
      const full=[{l:`Moles of ${fx(r0.f)}`,a:nR,u:'mol',w:`n = m ÷ Mr = ${fmt(mass)} ÷ ${fmt(M)} = ${fmt(nR)} mol`},{l:`Moles of ${fx(P.f)} (theoretical)`,a:nP,u:'mol',w:`Use the ratio ${R.coefs[0]} : ${R.coefs[nr+j]}. n = ${fmt(nR)} × ${R.coefs[nr+j]} ÷ ${R.coefs[0]} = ${fmt(nP)} mol`},{l:'Theoretical yield',a:theo,u,w:P.form==='gas'?`V = n × 24 = ${fmt(nP)} × 24 = ${fmt(theo)} dm³`:`m = n × Mr = ${fmt(nP)} × ${fmt(mr(P.f))} = ${fmt(theo)} g`},{l:'Percentage yield',a:pp,u:'%',w:`% yield = ${fmt(act)} ÷ ${fmt(theo)} × 100 = ${fmt(pp)} %`}];
      Y.q={eq,t:`${mass} g of ${r0.name} (${fx(r0.f)}) reacts completely. Any other reactant is in excess. The experiment collects ${fmt(act)} ${u} of ${P.name} (${fx(P.f)}). Calculate the theoretical yield and the percentage yield.`,b:S.mode==='challenge'?full.slice(2):full};
      break;
    }
  }
  Y.vals={};Y.checked=false;Y.shown=false;Y.hint=false;Y.done=false;
}
const yRight=(i)=>{const v=Y.vals[i];if(v===undefined||v==='')return null;const x=parseFloat(String(v).replace(',','.')),a=Y.q.b[i].a;return !isNaN(x)&&Math.abs(x-a)<=Math.abs(a)*0.02+1e-9;};
function renderYield(){
  if(!Y.q)newY();
  const m=S.mode,q=Y.q,w=m==='worked',rows=q.b.map((b,i)=>{
    const st=Y.shown||w?'shown':(Y.checked||m==='guided')&&yRight(i)!==null?(yRight(i)?'ok':'bad'):'';
    const val=Y.shown||w?fmt(b.a):(Y.vals[i]!==undefined?Y.vals[i]:'');
    return `<div class="frow"><div><b>${b.l}</b></div><div><input class="blank ${st}" data-yin="${i}" type="text" inputmode="decimal" autocomplete="off" value="${esc(val)}" aria-label="${esc(b.l.replace(/<[^>]*>/g,''))}"${w?' disabled':''}> <span class="hint">${b.u}</span></div><div class="hint">${(w||Y.shown)&&b.w?b.w:''}</div></div>`;
  }).join('');
  const formula=(w||m==='guided'||Y.hint)?`<div class="callout"><span class="tri" style="display:flex;flex-direction:column;gap:3px"><span>% yield = actual yield ÷ theoretical yield × 100</span><span>actual = theoretical × % yield ÷ 100</span><span>theoretical = actual ÷ % yield × 100</span></span></div>`:'';
  const simpleW=(w||Y.shown)&&Y.kind==='simple'?`<div class="callout"><b>Working.</b> ${q.w}</div>`:'';
  $('yView').innerHTML=`<section class="stage"><header><h2>Percentage yield practice</h2><p>The theoretical yield is the most product the equation allows. The actual yield is what you really collect. Percentage yield compares them. Answers within 2% are accepted.</p></header>
  <div class="chips" role="group" aria-label="Question type"><button class="chip" data-act="ykind" data-k="simple" aria-pressed="${Y.kind==='simple'}">Quick calculations</button><button class="chip" data-act="ykind" data-k="eq" aria-pressed="${Y.kind==='eq'}">With an equation</button></div>
  <div class="q">${q.eq?`<div class="eq chem"><span>${q.eq}</span></div>`:''}<div class="text">${q.t}</div>${formula}<div class="rows">${rows}</div>${simpleW}
   <div class="chips">${w?'':'<button class="btn teal" data-act="ycheck">Check</button>'}${m==='standard'&&!Y.hint&&!w?'<button class="btn" data-act="yhint">Hint</button>':''}${w?'':'<button class="btn" data-act="yshow">Show working</button>'}<button class="btn pri" data-act="ynew">New question</button></div>
   <div id="yres"></div></div>
  <div class="score">${Y.ok} correct out of ${Y.n} checks</div></section>`;
}

/* ---------- events ---------- */
function applyShow(ids){
  ids.forEach(id=>{
    if(!S.fields[id])return;
    S.shown[id]=true;S.vals[id]=ansStr(id);
    const el=document.querySelector(`[data-f="${id}"]`);if(el)el.value=S.vals[id];
  });
  paint();
}
const idsBy=p=>Object.keys(S.fields).filter(id=>p==='f'?/^f\d/.test(id):id.indexOf(p)===0);
document.addEventListener('click',e=>{
  const t=e.target.closest('[data-act]');if(!t)return;
  const d=t.dataset,a=d.act;
  if(a==='mode')setMode(d.m);
  else if(a==='tab'){S.tab=d.t;renderAll();}
  else if(a==='cat'){S.cat=d.c;setReaction();}
  else if(a==='random'){
    for(let k=0;k<40;k++){const p={};if(S.cat==='neutralisation'){p.type=rnd(Object.keys(NTYPE))[0]?rnd(Object.keys(NTYPE)):'hyd';p.acid=rnd(OPTS.neutralisation[0].list)[0];p.base=rnd(NBASE[p.type]);}else OPTS[S.cat].forEach(o=>{p[o.k]=rnd(o.list)[0];});S.pick[S.cat]=p;if(build(S.cat,p).valid)break;}
    setReaction();
  }
  else if(a==='go')go(+d.s);
  else if(a==='show')applyShow([d.id]);
  else if(a==='showall')applyShow(idsBy(d.p));
  else if(a==='showrow'){const rm=S.rowMeta.find(r=>r.k===d.r);if(rm)applyShow(rm.ids);}
  else if(a==='check'){idsBy(d.p).forEach(id=>{S.checked[id]=true;});paint();}
  else if(a==='checkrow'){const rm=S.rowMeta.find(r=>r.k===d.r);if(rm){rm.ids.forEach(id=>{S.checked[id]=true;});paint();}}
  else if(a==='hint'){S.hintOpen[d.id]=!S.hintOpen[d.id];if(d.id==='atoms')refreshBal();else renderStage();}
  else if(a==='coef'){const i=+d.i;S.coefs[i]=Math.max(1,Math.min(30,S.coefs[i]+(+d.d)));$('co'+i).textContent=S.coefs[i];S.checked.bal=false;refreshBal();}
  else if(a==='checkbal'){S.checked.bal=true;refreshBal();}
  else if(a==='solvebal'){S.coefs=S.R.coefs.slice();S.checked.bal=true;renderEq();renderRail();renderStage();}
  else if(a==='adj'){adj(+d.i,d.k,+d.d);}
  else if(a==='randamt'){S.A=defaults(S.R);if(S.prac==='tit'){S.prac='mix';S.tit=null;}clearCalc();renderStage();}
  else if(a==='run'){const c=calc(S.R,S.A);S.run={c};$('runout').innerHTML=runHTML();renderRail();}
  else if(a==='prac'){S.prac=d.p;if(d.p==='tit')setupTit();else{S.tit=null;S.A=defaults(S.R);clearCalc();}renderRail();renderStage();}
  else if(a==='titadd')titAdd(+d.d);
  else if(a==='titreset'){S.tit.added=0;S.tit.msg='';S.tit.done=false;titRefresh();}
  else if(a==='titrecord')titRecord();
  else if(a==='revnext'){S.reveal++;renderStage();}
  else if(a==='revall'){S.reveal=99;renderStage();}
  else if(a==='revreset'){S.reveal=0;renderStage();}
  else if(a==='bfilter'){B.filter=d.k;newBal();renderBal();}
  else if(a==='blevel'){B.level=d.k;newBal();renderBal();}
  else if(a==='bnew'){newBal();renderBal();}
  else if(a==='bco'){const i=+d.i;B.co[i]=Math.max(1,Math.min(30,B.co[i]+(+d.d)));B.checked=false;renderBal();}
  else if(a==='bcheck'){B.checked=true;B.n++;if(balInfoB().every(x=>x.ok)&&B.co.reduce(gcd)===1&&!B.done&&!B.shown){B.ok++;B.done=true;}renderBal();}
  else if(a==='bhint'){B.hint=true;renderBal();}
  else if(a==='bshow'){B.co=B.R.coefs.slice();B.shown=true;renderBal();}
  else if(a==='ykind'){Y.kind=d.k;Y.ok=0;Y.n=0;newY();renderYield();}
  else if(a==='ynew'){newY();renderYield();}
  else if(a==='yhint'){Y.hint=true;renderYield();}
  else if(a==='yshow'){Y.shown=true;renderYield();}
  else if(a==='ycheck'){
    const all=Y.q.b.every((b,i)=>yRight(i)===true);
    if(Y.q.b.some((b,i)=>yRight(i)===null)){Y.checked=true;renderYield();$('yres').innerHTML='<div class="callout warn">Fill in every blank first.</div>';}
    else{Y.checked=true;Y.n++;if(all&&!Y.done&&!Y.shown){Y.ok++;Y.done=true;}renderYield();$('yres').innerHTML=all?'<div class="callout good"><b>Correct.</b></div>':'<div class="callout bad">Not quite. Check each step, your units and your rounding.</div>';}
  }
  else if(a==='tkind'){T.kind=d.k;newQ();T.ok=0;T.n=0;renderTool();}
  else if(a==='tnew'){newQ();renderTool();}
  else if(a==='tcheck'){toolCheck();$('tres').innerHTML=toolRes();const sc=document.querySelector('#toolView .score');if(sc)sc.textContent=`${T.ok} correct out of ${T.n} attempts`;}
  else if(a==='thint'){T.hint=true;renderTool();}
  else if(a==='tshow'){T.shown=true;$('tres').innerHTML=toolRes();}
});
function onField(el){
  const id=el.dataset.f;S.vals[id]=el.value;S.shown[id]=false;
  if(S.mode!=='guided')S.checked[id]=false;
  paint();
}
document.addEventListener('input',e=>{
  const el=e.target;
  if(el.dataset.f!==undefined&&el.tagName==='INPUT')onField(el);
  else if(el.dataset.amt!==undefined){const i=+el.dataset.amt,v=parseFloat(el.value);if(!isNaN(v)&&v>=0){S.A[i][el.dataset.k]=v;$('vis'+i).innerHTML=visHTML(i);clearCalc();}}
  else if(el.dataset.tin!==undefined){T.val=el.value;}
  else if(el.dataset.yin!==undefined){Y.vals[+el.dataset.yin]=el.value;Y.checked=false;el.classList.remove('ok','bad');}
});
document.addEventListener('change',e=>{
  const el=e.target;
  if(el.dataset.f!==undefined&&el.tagName==='SELECT')onField(el);
  else if(el.dataset.pick){const pk=S.pick[S.cat];pk[el.dataset.pick]=el.value;if(S.cat==='neutralisation'&&el.dataset.pick==='type')pk.base=NBASE[el.value][0];setReaction();}
  else if(el.dataset.amtC!==undefined){const i=+el.dataset.amtC;S.A[i].c=parseFloat(el.value);$('vis'+i).innerHTML=visHTML(i);clearCalc();}
  else if(el.dataset.titInd!==undefined){S.tit.ind=el.value;titRefresh();}
});
let holdT=null;
document.addEventListener('pointerdown',e=>{
  const b=e.target.closest('[data-hold]');if(!b)return;
  const i=+b.dataset.i,k=b.dataset.k,d=+b.dataset.d;
  adj(i,k,d);clearInterval(holdT);holdT=setInterval(()=>adj(i,k,d),90);
});
['pointerup','pointercancel','pointerleave','blur'].forEach(ev=>document.addEventListener(ev,()=>{clearInterval(holdT);holdT=null;},true));
document.addEventListener('keydown',e=>{
  if((e.key==='Enter'||e.key===' ')&&e.target.matches&&e.target.matches('[data-hold]')){e.preventDefault();adj(+e.target.dataset.i,e.target.dataset.k,+e.target.dataset.d);}
  if(e.key==='Enter'&&e.target.id==='tin'){toolCheck();$('tres').innerHTML=toolRes();const sc=document.querySelector('#toolView .score');if(sc)sc.textContent=`${T.ok} correct out of ${T.n} attempts`;}
});
/* ---------- floating calculator ---------- */
(function(){
  const box=$('calc'),inp=$('calcIn'),hist=$('calcHist'),btn=$('calcBtn'),bar=$('calcBar');
  let ans=0,pos=null;
  const store={get(k){try{return localStorage.getItem(k);}catch(e){return null;}},set(k,v){try{localStorage.setItem(k,v);}catch(e){}}};
  function evaluate(src){
    const s=src.replace(/\s/g,'').replace(/×/g,'*').replace(/÷/g,'/').replace(/−/g,'-').replace(/Ans/g,'#');
    const re=/\d+\.?\d*(?:e[+-]?\d+)?|\.\d+(?:e[+-]?\d+)?|[()+\-*\/^]|√|#/g;
    const toks=s.match(re)||[];
    if(toks.join('')!==s)throw new Error('bad');
    let p=0;
    const peek=()=>toks[p],next=()=>toks[p++];
    function expr(){let v=term();while(peek()==='+'||peek()==='-'){const o=next(),r=term();v=o==='+'?v+r:v-r;}return v;}
    function term(){let v=pow();for(;;){const t=peek();if(t==='*'||t==='/'){next();const r=pow();v=t==='*'?v*r:v/r;}else if(t==='('||t==='√'||t==='#'||(t&&/^[\d.]/.test(t))){v*=pow();}else break;}return v;}
    function pow(){const b=unary();if(peek()==='^'){next();return Math.pow(b,pow());}return b;}
    function unary(){if(peek()==='-'){next();return -unary();}if(peek()==='+'){next();return unary();}return atom();}
    function atom(){const t=next();if(t===undefined)throw new Error('bad');if(t==='('){const v=expr();if(next()!==')')throw new Error('bad');return v;}if(t==='√'){const v=atom();if(v<0)throw new Error('bad');return Math.sqrt(v);}if(t==='#')return ans;const n=parseFloat(t);if(isNaN(n))throw new Error('bad');return n;}
    const v=expr();if(p!==toks.length)throw new Error('bad');
    if(!isFinite(v))throw new Error('bad');
    return v;
  }
  const show=v=>{const a=Math.abs(v);return String((a!==0&&(a<1e-4||a>=1e10))?v.toExponential(6).replace(/\.?0+e/,'e'):parseFloat(v.toPrecision(10)));};
  function equals(){
    const src=inp.value;if(!src.trim())return;
    try{const v=evaluate(src);ans=v;hist.textContent=src+' =';inp.value=show(v);}
    catch(e){hist.textContent='Check the brackets and symbols';}
  }
  function key(k){
    if(k==='C'){inp.value='';hist.textContent='';}
    else if(k==='BS')inp.value=inp.value.replace(/(Ans|.)$/,'');
    else if(k==='=')equals();
    else{
      const st=inp.selectionStart,en=inp.selectionEnd,v=inp.value;
      if(document.activeElement===inp&&st!=null){inp.value=v.slice(0,st)+k+v.slice(en);const c=st+k.length;inp.setSelectionRange(c,c);}
      else inp.value=v+k;
    }
  }
  function clamp(){
    if(!pos)return;
    const w=box.offsetWidth,h=box.offsetHeight;
    pos.x=Math.max(8,Math.min(window.innerWidth-w-8,pos.x));pos.y=Math.max(8,Math.min(window.innerHeight-h-8,pos.y));
    box.style.left=pos.x+'px';box.style.top=pos.y+'px';
  }
  function setOn(on){
    box.hidden=!on;btn.setAttribute('aria-pressed',on);btn.textContent='Calculator: '+(on?'on':'off');
    btn.className=on?'btn teal':'btn';store.set('moleLabCalc',on?'1':'0');
    if(on){
      if(!pos){let sx=null;try{sx=JSON.parse(store.get('moleLabCalcPos'));}catch(e){}
        pos=sx&&typeof sx.x==='number'?sx:{x:window.innerWidth-box.offsetWidth-16,y:window.innerHeight-box.offsetHeight-16};}
      clamp();
    }
  }
  function setType(sci){
    box.classList.toggle('sci',sci);
    box.querySelectorAll('[data-cm]').forEach(b=>b.setAttribute('aria-pressed',(b.dataset.cm==='sci')===sci));
    store.set('moleLabCalcType',sci?'sci':'basic');clamp();
  }
  btn.addEventListener('click',()=>{if(box.hidden&&window.innerWidth<700&&!$('ptp').hidden)$('ptClose').click();setOn(box.hidden);});
  $('calcClose').addEventListener('click',()=>setOn(false));
  $('calcMin').addEventListener('click',()=>{box.classList.toggle('min');$('calcMin').textContent=box.classList.contains('min')?'Open':'Min';clamp();});
  $('calcKeys').addEventListener('click',e=>{const b=e.target.closest('[data-k]');if(b)key(b.dataset.k);});
  box.querySelectorAll('[data-cm]').forEach(b=>b.addEventListener('click',()=>setType(b.dataset.cm==='sci')));
  inp.addEventListener('keydown',e=>{if(e.key==='Enter'){e.preventDefault();equals();}else if(e.key==='Escape'){setOn(false);btn.focus();}});
  inp.addEventListener('input',()=>{inp.value=inp.value.replace(/\*/g,'×').replace(/\//g,'÷').replace(/-/g,'−');});
  let drag=null;
  bar.addEventListener('pointerdown',e=>{
    if(e.target.closest('button'))return;
    drag={dx:e.clientX-pos.x,dy:e.clientY-pos.y};bar.setPointerCapture(e.pointerId);
  });
  bar.addEventListener('pointermove',e=>{if(!drag)return;pos.x=e.clientX-drag.dx;pos.y=e.clientY-drag.dy;clamp();});
  const end=()=>{if(drag){drag=null;store.set('moleLabCalcPos',JSON.stringify(pos));}};
  bar.addEventListener('pointerup',end);bar.addEventListener('pointercancel',end);
  window.addEventListener('resize',clamp);
  if(store.get('moleLabCalcType')==='sci')setType(true);
  if(store.get('moleLabCalc')==='1')setOn(true);
})();

/* ---------- periodic table panel ---------- */
(function(){
  const EL=`H,Hydrogen,1;He,Helium,4;Li,Lithium,7;Be,Beryllium,9;B,Boron,11;C,Carbon,12;N,Nitrogen,14;O,Oxygen,16;F,Fluorine,19;Ne,Neon,20;Na,Sodium,23;Mg,Magnesium,24;Al,Aluminium,27;Si,Silicon,28;P,Phosphorus,31;S,Sulfur,32;Cl,Chlorine,35.5;Ar,Argon,40;K,Potassium,39;Ca,Calcium,40;Sc,Scandium,45;Ti,Titanium,48;V,Vanadium,51;Cr,Chromium,52;Mn,Manganese,55;Fe,Iron,56;Co,Cobalt,59;Ni,Nickel,59;Cu,Copper,64;Zn,Zinc,65;Ga,Gallium,70;Ge,Germanium,73;As,Arsenic,75;Se,Selenium,79;Br,Bromine,80;Kr,Krypton,84;Rb,Rubidium,85;Sr,Strontium,88;Y,Yttrium,89;Zr,Zirconium,91;Nb,Niobium,93;Mo,Molybdenum,96;Tc,Technetium,98;Ru,Ruthenium,101;Rh,Rhodium,103;Pd,Palladium,106;Ag,Silver,108;Cd,Cadmium,112;In,Indium,115;Sn,Tin,119;Sb,Antimony,122;Te,Tellurium,128;I,Iodine,127;Xe,Xenon,131;Cs,Caesium,133;Ba,Barium,137;La,Lanthanum,139;Ce,Cerium,140;Pr,Praseodymium,141;Nd,Neodymium,144;Pm,Promethium,147;Sm,Samarium,150;Eu,Europium,152;Gd,Gadolinium,157;Tb,Terbium,159;Dy,Dysprosium,163;Ho,Holmium,165;Er,Erbium,167;Tm,Thulium,169;Yb,Ytterbium,173;Lu,Lutetium,175;Hf,Hafnium,178;Ta,Tantalum,181;W,Tungsten,184;Re,Rhenium,186;Os,Osmium,190;Ir,Iridium,192;Pt,Platinum,195;Au,Gold,197;Hg,Mercury,201;Tl,Thallium,204;Pb,Lead,207;Bi,Bismuth,209;Po,Polonium,209;At,Astatine,210;Rn,Radon,222;Fr,Francium,223;Ra,Radium,226;Ac,Actinium,227;Th,Thorium,232;Pa,Protactinium,231;U,Uranium,238;Np,Neptunium,237;Pu,Plutonium,244;Am,Americium,243;Cm,Curium,247;Bk,Berkelium,247;Cf,Californium,251;Es,Einsteinium,252;Fm,Fermium,257;Md,Mendelevium,258;No,Nobelium,259;Lr,Lawrencium,262`.split(';').map((s,i)=>{const [sym,name,ar]=s.split(',');return {z:i+1,sym,name,ar:parseFloat(ar)};});
  const box=$('ptp'),btn=$('ptBtn'),bar=$('ptBar'),view=$('ptView'),det=$('ptDetail'),search=$('ptSearch');
  const store={get(k){try{return localStorage.getItem(k);}catch(e){return null;}},set(k,v){try{localStorage.setItem(k,v);}catch(e){}}};
  const DEF='Tap an element to pin its atomic number, symbol and relative atomic mass here. It stays until you choose another.';
  const picks=[];let multi=false;
  const has=z=>picks.some(p=>p.z===z),cnt=z=>{const p=picks.find(p=>p.z===z);return p?p.n:0;};
  const SYM={};EL.forEach(e=>{SYM[e.sym]=e;});
  function addPick(z,d){
    let p=picks.find(x=>x.z===z);
    if(!p){if(d<=0)return;p={z,n:0};picks.push(p);}
    p.n=Math.min(99,p.n+d);
    if(p.n<=0)picks.splice(picks.indexOf(p),1);
  }
  const fmtn=x=>String(parseFloat(Number(x).toPrecision(8)));
  function parseFormula(raw){
    let s=raw.replace(/[₀-₉]/g,d=>'₀₁₂₃₄₅₆₇₈₉'.indexOf(d)).replace(/\s+/g,'').replace(/\[/g,'(').replace(/\]/g,')');
    if(!s)return null;
    const total={};
    for(const seg of s.split(/[.·•*]/)){
      const mm=seg.match(/^(\d*)(.*)$/),mult=mm[1]?parseInt(mm[1],10):1,body=mm[2];
      if(!body)throw new Error('Something is missing after the dot');
      let i=0;
      const num=()=>{let t='';while(i<body.length&&/\d/.test(body[i]))t+=body[i++];return t?parseInt(t,10):1;};
      const grp=()=>{
        const m={};
        while(i<body.length){
          const c=body[i];
          if(c==='('){i++;const sub=grp();if(body[i]!==')')throw new Error('A bracket is not closed');i++;const n=num();for(const k in sub)m[k]=(m[k]||0)+sub[k]*n;}
          else if(c===')')return m;
          else if(/[A-Z]/.test(c)){let el=c;i++;while(i<body.length&&/[a-z]/.test(body[i]))el+=body[i++];if(!SYM[el])throw new Error('“'+el+'” is not an element symbol');m[el]=(m[el]||0)+num();}
          else if(/[a-z]/.test(c))throw new Error('Element symbols start with a capital letter, for example Cu, not cu');
          else throw new Error('Unexpected “'+c+'”');
        }
        return m;
      };
      const r=grp();
      if(i<body.length)throw new Error('There is a ) without a matching (');
      for(const k in r)total[k]=(total[k]||0)+r[k]*mult;
    }
    return total;
  }
  function renderMr(){
    const out=$('ptMrOut'),raw=$('ptFormula').value;
    if(!raw.trim()){out.innerHTML='';return;}
    try{
      const t=parseFormula(raw),keys=Object.keys(t);
      if(!keys.length){out.innerHTML='';return;}
      let sum=0;
      const rows=keys.map(k=>{const sub=t[k]*SYM[k].ar;sum+=sub;return `<tr><td><b>${k}</b></td><td>${t[k]} × ${SYM[k].ar}</td><td>= ${fmtn(sub)}</td></tr>`;}).join('');
      out.innerHTML=`<table><tbody>${rows}</tbody></table><div class="res">Mr = <b>${fmtn(sum)}</b></div>`;
    }catch(err){out.innerHTML=`<div class="callout warn" style="padding:6px 10px">${err.message}</div>`;}
  }
  let pos=null,mode='table',sel=null,size=store.get('moleLabPTSize')||(window.innerWidth<900?'s':'m'),fade=false;
  function applySize(){
    ['s','m','l'].forEach(k=>box.classList.toggle('sz-'+k,k===size));
    box.querySelectorAll('[data-sz]').forEach(b=>b.setAttribute('aria-pressed',b.dataset.sz===size));
    store.set('moleLabPTSize',size);
  }
  function applyFade(){box.classList.toggle('fade',fade);$('ptFade').setAttribute('aria-pressed',fade);$('ptFade').textContent=fade?'Fade: on, tap to restore':'Fade';}
  const inList=(z,a)=>a.includes(z);
  function cat(z){
    if(inList(z,[3,11,19,37,55,87]))return ['k-alk','Group 1 metal (alkali metal)'];
    if(inList(z,[4,12,20,38,56,88]))return ['k-aea','Group 2 metal (alkaline earth metal)'];
    if(inList(z,[2,10,18,36,54,86]))return ['k-nob','Noble gas (Group 18)'];
    if(inList(z,[9,17,35,53,85]))return ['k-hal','Halogen (Group 17)'];
    if(z>=57&&z<=71)return ['k-lan','Lanthanoid'];
    if(z>=89)return ['k-act','Actinoid'];
    if((z>=21&&z<=30)||(z>=39&&z<=48)||(z>=72&&z<=80))return ['k-trn','Transition metal'];
    if(inList(z,[13,31,49,50,81,82,83]))return ['k-pt','Other metal'];
    if(inList(z,[5,14,32,33,51,52,84]))return ['k-mtl','Metalloid'];
    return ['k-nm','Non-metal'];
  }
  function place(z){
    if(z===1)return [1,2];if(z===2)return [18,2];
    if(z<=4)return [z-2,3];if(z<=10)return [13+z-5,3];
    if(z<=12)return [z-10,4];if(z<=18)return [13+z-13,4];
    if(z<=36)return [z-18,5];if(z<=54)return [z-36,6];
    if(z<=56)return [z-54,7];if(z<=71)return [3+z-57,10];
    if(z<=86)return [4+z-72,7];if(z<=88)return [z-86,8];
    return [3+z-89,11];
  }
  function info(e){
    const [c,r]=place(e.z),k=cat(e.z);
    const gp=(e.z>=57&&e.z<=71)||(e.z>=89)?'':`Group ${c}, `;
    const per=r===10?6:r===11?7:r-1;
    return `<b>${e.sym}</b> ${e.name} · atomic number <b>${e.z}</b> · relative atomic mass (Ar) <b>${e.ar}</b><br><span style="color:var(--muted);font-size:13px">${gp}Period ${per} · ${k[1]}</span>`;
  }
  function match(e){
    const q=search.value.trim().toLowerCase();
    if(!q)return true;
    return e.sym.toLowerCase()===q||e.name.toLowerCase().includes(q)||String(e.z)===q||e.sym.toLowerCase().startsWith(q);
  }
  function renderPicks(){
    const pk=$('ptPicks');
    if(!picks.length&&!multi){pk.hidden=true;pk.innerHTML='';return;}
    pk.hidden=false;
    if(!picks.length){pk.innerHTML='<div>Select several is on. Tap elements to add them. Tap the same element again to add another of it.</div>';return;}
    let sum=0;const terms=picks.map(p=>{const e=EL[p.z-1];sum+=p.n*e.ar;return (p.n>1?p.n+' × ':'')+e.ar;});
    const total=picks.reduce((t,p)=>t+p.n,0);
    pk.innerHTML=`<div><b>My list</b> (${total} atom${total===1?'':'s'})</div><div class="row">${picks.map(p=>{const e=EL[p.z-1];return `<span class="pk"><b>${e.sym}</b> ×${p.n} · Ar ${e.ar}<button data-dec="${p.z}" aria-label="One fewer ${e.name}">−</button><button data-inc="${p.z}" aria-label="One more ${e.name}">+</button><button class="x" data-rm="${p.z}" aria-label="Remove ${e.name} from the list">✕</button></span>`;}).join('')}<button class="chip" data-clearpicks style="background:var(--bad-soft);color:var(--bad)">Clear list</button></div><div class="tot">Total: ${terms.join(' + ')} = <b>${fmtn(sum)}</b>. This is the Mr if your list is the whole formula.</div>`;
  }
  function render(){
    const q=search.value.trim();
    if(mode==='list'){
      const rows=EL.filter(match).map(e=>`<tr data-z="${e.z}" style="cursor:pointer${has(e.z)?';background:var(--teal-soft)':sel===e.z?';background:var(--accent-soft)':''}"><td>${e.z}</td><td><b>${e.sym}</b></td><td>${e.name}</td><td>${e.ar}</td></tr>`).join('');
      view.innerHTML=`<div class="ptl"><table><thead><tr><th>Atomic number</th><th>Symbol</th><th>Name</th><th>Ar</th></tr></thead><tbody>${rows||'<tr><td colspan="4">No element matches.</td></tr>'}</tbody></table></div>`;
      return;
    }
    let h='<div class="pt" role="group" aria-label="Periodic table of elements">';
    for(let g=1;g<=18;g++)h+=`<div class="g" style="grid-column:${g};grid-row:1">${g}</div>`;
    h+=`<div class="ph" style="grid-column:3;grid-row:7">57–71</div><div class="ph" style="grid-column:3;grid-row:8">89–103</div>`;
    EL.forEach(e=>{
      const [c,r]=place(e.z),k=cat(e.z)[0],dim=q&&!match(e);
      h+=`<button class="el ${k}${dim?' dim':''}${sel===e.z?' sel':''}${has(e.z)?' pick':''}" data-z="${e.z}" style="grid-column:${c};grid-row:${r}" aria-label="${e.name}, atomic number ${e.z}, relative atomic mass ${e.ar}"><span class="z">${e.z}</span><span class="s">${e.sym}</span><span class="a">${e.ar}</span>${cnt(e.z)>1?`<span class="cnt">×${cnt(e.z)}</span>`:''}</button>`;
    });
    view.innerHTML=h+'</div>';
  }
  function showDetail(z){sel=z;det.innerHTML=info(EL[z-1])+' <button class="chip" data-clr style="padding:2px 10px;font-size:12px;margin-left:6px">Clear</button>';view.querySelectorAll('.el.sel').forEach(b=>b.classList.remove('sel'));const b=view.querySelector(`[data-z="${z}"]`);if(b)b.classList.add('sel');}
  view.addEventListener('click',e=>{
    const b=e.target.closest('[data-z]');if(!b)return;
    const z=+b.dataset.z;
    if(multi){addPick(z,1);render();renderPicks();showDetail(z);}
    else showDetail(z);
  });
  $('ptPicks').addEventListener('click',e=>{
    const inc=e.target.closest('[data-inc]'),dec=e.target.closest('[data-dec]'),rm=e.target.closest('[data-rm]');
    if(inc)addPick(+inc.dataset.inc,1);
    else if(dec)addPick(+dec.dataset.dec,-1);
    else if(rm){const p=picks.find(x=>x.z===+rm.dataset.rm);if(p)picks.splice(picks.indexOf(p),1);}
    else if(e.target.closest('[data-clearpicks]'))picks.length=0;
    else return;
    render();renderPicks();
  });
  $('ptFormula').addEventListener('input',renderMr);
  $('ptMulti').addEventListener('click',()=>{
    multi=!multi;$('ptMulti').setAttribute('aria-pressed',multi);$('ptMulti').textContent=multi?'Select several: on':'Select several';
    renderPicks();
  });
  det.addEventListener('click',e=>{if(e.target.closest('[data-clr]')){sel=null;view.querySelectorAll('.el.sel').forEach(x=>x.classList.remove('sel'));det.textContent=DEF;}});
  view.addEventListener('mouseover',e=>{if(sel!==null)return;const b=e.target.closest('[data-z]');if(b)det.innerHTML=info(EL[+b.dataset.z-1]);});
  view.addEventListener('focusin',e=>{if(sel!==null)return;const b=e.target.closest('[data-z]');if(b)det.innerHTML=info(EL[+b.dataset.z-1]);});
  search.addEventListener('input',()=>{
    render();
    const m=EL.filter(match);
    if(search.value.trim()&&m.length===1){showDetail(m[0].z);}
    else if(!search.value.trim()&&sel===null){det.textContent=DEF;}
  });
  box.querySelectorAll('[data-sz]').forEach(b=>b.addEventListener('click',()=>{size=b.dataset.sz;applySize();clamp();}));
  $('ptFade').addEventListener('click',()=>{fade=!fade;applyFade();});
  box.querySelectorAll('[data-pv]').forEach(b=>b.addEventListener('click',()=>{
    mode=b.dataset.pv;box.querySelectorAll('[data-pv]').forEach(x=>x.setAttribute('aria-pressed',x.dataset.pv===mode));render();clamp();
  }));
  function clamp(){
    if(!pos)return;
    const w=box.offsetWidth,hh=box.offsetHeight;
    pos.x=Math.max(8,Math.min(window.innerWidth-w-8,pos.x));pos.y=Math.max(8,Math.min(window.innerHeight-hh-8,pos.y));
    box.style.left=pos.x+'px';box.style.top=pos.y+'px';
  }
  function setOn(on){
    box.hidden=!on;btn.setAttribute('aria-pressed',on);btn.textContent='Periodic table: '+(on?'on':'off');
    btn.className=on?'btn teal':'btn';store.set('moleLabPT',on?'1':'0');
    if(on){
      if(!view.innerHTML)render();renderPicks();
      if(!pos){let sx=null;try{sx=JSON.parse(store.get('moleLabPTPos'));}catch(e){}
        pos=sx&&typeof sx.x==='number'?sx:{x:16,y:90};}
      clamp();
    }
  }
  btn.addEventListener('click',()=>{if(box.hidden&&window.innerWidth<700&&!$('calc').hidden)$('calcClose').click();setOn(box.hidden);});
  $('ptClose').addEventListener('click',()=>setOn(false));
  $('ptMin').addEventListener('click',()=>{box.classList.toggle('min');$('ptMin').textContent=box.classList.contains('min')?'Full':'Compact';if(box.classList.contains('min')){fade=false;applyFade();}clamp();});
  let drag=null;
  bar.addEventListener('pointerdown',e=>{if(e.target.closest('button'))return;drag={dx:e.clientX-pos.x,dy:e.clientY-pos.y};bar.setPointerCapture(e.pointerId);});
  bar.addEventListener('pointermove',e=>{if(!drag)return;pos.x=e.clientX-drag.dx;pos.y=e.clientY-drag.dy;clamp();});
  const end=()=>{if(drag){drag=null;store.set('moleLabPTPos',JSON.stringify(pos));}};
  bar.addEventListener('pointerup',end);bar.addEventListener('pointercancel',end);
  window.addEventListener('resize',clamp);
  applySize();
  if(store.get('moleLabPT')==='1')setOn(true);
})();

setReaction();
