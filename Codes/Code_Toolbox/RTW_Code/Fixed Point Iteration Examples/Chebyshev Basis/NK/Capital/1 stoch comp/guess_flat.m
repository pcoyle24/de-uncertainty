function pf = guess_flat(O,P,S,G)

pf.pi = ones(G.nl.griddim)*P.pi;
pf.n = ones(G.nl.griddim)*P.n;
pf.k = ones(G.nl.griddim)*S.k;
end