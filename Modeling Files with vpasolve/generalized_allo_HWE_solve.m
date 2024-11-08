

syms g00 g01 g10 g11 s h1 h2 h3 mu nu beta gamma

assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);

sel_g00 = 0 == g00^2 + (.5+(beta+beta*gamma)/32)*2*g01*g00 + (.5+(beta+beta*gamma)/32)*2*g10*g00 + ((3*beta+beta*gamma)/16)*g01^2 + ((3*beta+beta*gamma)/16)*g10^2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + ((beta+beta*gamma)/32)*2*g01*g11 + ((beta+beta*gamma)/32)*2*g10*g11;
sel_g10 = 0 == ((3*beta+3*beta*gamma)/32)*2*g01*g00 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g00 + ((beta+3*beta*gamma)/16)*g01^2 + (1-(7*beta+5*beta*gamma)/16)*g10^2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + ((3*beta+3*beta*gamma)/32)*2*g01*g11 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g11;
sel_g01 = 0 == (.5-(5*beta+5*beta*gamma)/32)*2*g01*g00 + ((3*beta+3*beta*gamma)/32)*2*g10*g00 + (1-(7*beta+5*beta*gamma)/16)*g01^2 + ((beta+3*beta*gamma)/16)*g10^2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + (.5-(5*beta+5*beta*gamma)/32)*2*g01*g11 + ((3*beta+3*beta*gamma)/32)*2*g10*g11;
sel_g11 = 0 == ((beta+beta*gamma)/32)*2*g01*g00 + ((beta+beta*gamma)/32)*2*g10*g00 + ((3*beta+beta*gamma)/16)*g01^2 + ((3*beta+beta*gamma)/16)*g10^2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + (.5+(beta+beta*gamma)/32)*2*g01*g11 + (.5+(beta+beta*gamma)/32)*2*g10*g11 + g11^2;

eqn5 = g00+g01+g10+g11 == 1;

S = solve(sel_g01, sel_g10, sel_g11, eqn5, g00, g01, g10, g11, 'ReturnConditions', true);