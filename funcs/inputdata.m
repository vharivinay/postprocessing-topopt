function [nelx,nely,volfrac,penal,rmin] = inputdata(specimen,sf)

if nargin < 2
    sf = 1.0;
end

switch specimen
    case 'MBB'
        nelx = 150*sf;
        nely = 50*sf;
        volfrac=0.5;
        penal=3;
        rmin=0.03*nelx;
    case 'half-wheel'
        nelx = 120*sf;
        nely = 50*sf;
        volfrac=0.5;
        penal=3;
        rmin=0.03*nelx;
    case {'short_cantilever_endload','short_cantilever_midload'}
        nelx = 160*sf;
        nely = 100*sf;
        volfrac=0.4;
        penal=3;
        rmin=0.03*nelx;
    case {'plate_w_hole_endload','plate_w_hole_midload'}
        nelx = 160*sf;
        nely = 100*sf;
        volfrac=0.4;
        penal=3;
        rmin=0.03*nelx;
    case 'multiload_cantilever'
        nelx = 100*sf;
        nely = 100*sf;
        volfrac=0.4;
        penal=3;
        rmin=0.03*nelx;
    case 'l-beam'
        nelx = 100*sf;
        nely = 100*sf;
        volfrac=0.25;
        penal=3;
        rmin=0.03*nelx;
    case {'t-beam','t-beam_w_hole'}
        nelx = 120*sf;
        nely = 100*sf;
        volfrac=0.25;
        penal=3;
        rmin=0.03*nelx;
    otherwise
        error('inputdata:specimenNotFound',"\nAllowed Inputs are:\n'MBB','half-wheel','short_cantilever_endload','short_cantilever_midload','multiload_cantilever'\n'plate_w_hole_endlaod','plate_w_hole_midlaod','l-beam','t-beam',\n't-beam_w_hole'")
end