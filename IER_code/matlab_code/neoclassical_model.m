function [y, x, fx,fxp, fy,fyp, f] = neoclassical_model

syms KAPA ALFA_PAI ALFA ALFA_R ALFA_Y Y PAI TH BETTA R RHO_th RHO_eps RHO_gam EPS GAM 
syms y_cu y_cup pai_cu pai_cup r_ba1 r_ba1p pai_ba1 pai_ba1p r_cu r_cup
syms th_ba1  th_ba1p
syms eps_ba1 eps_ba1p
syms gam_ba1 gam_ba1p

f = [];
f = [f;- 1/y_cu          + BETTA*(1/y_cup) * (r_ba1p/pai_cup) + log(eps_ba1/EPS)];
f = [f;- log(pai_cu/PAI) + BETTA*log(pai_cup/PAI) + KAPA*log(y_cu/Y) + log(gam_ba1/GAM)];
f = [f;- log(r_ba1p/R)   + ALFA_R*log(r_ba1/R) + ALFA_PAI*log(pai_cu/PAI) + ALFA_Y*log(y_cu/Y) + log(th_ba1/TH)];

f = [f;log( th_ba1p /TH  ) - RHO_th  * log( th_ba1  /TH )];
f = [f;log( eps_ba1p/EPS ) - RHO_eps * log( eps_ba1 /EPS)];
f = [f;log( gam_ba1p/GAM ) - RHO_gam * log( gam_ba1 /GAM)];
f = [f;- r_ba1p + r_cu];

% controls
y      = [ r_cu  y_cu  pai_cu  ];
yp     = [ r_cup y_cup pai_cup ];
Y      = [ R     Y      PAI    ];
ny     = length(y) - nz;

% Shocks
x      = [ r_ba1  th_ba1  eps_ba1  gam_ba1  ];
xp     = [ r_ba1p th_ba1p eps_ba1p gam_ba1p ];
XX     = [ R      TH       EPS      GAM     ];
nx     = length(x);

% Create function f
f = [ f1; f2; f3; f4; f5; f6; f7];

vs     = [ x, y, xp, yp ];
VS     = [ XX Y XX Y];
f      = subs(f, vs, exp(vs));
fx     = jacobian(f,x);
fxp    = jacobian(f,xp);
fy     = jacobian(f,y);
fyp    = jacobian(f,yp);
n      = length(f);

f      = subs(f,   vs, log(VS));
fx     = subs(fx,  vs, log(VS));
fxp    = subs(fxp, vs, log(VS));
fy     = subs(fy,  vs, log(VS));
fyp    = subs(fyp, vs, log(VS));

% find non-zero elements of f, fx, fy,fxp and fyp to exclude from evaluation
f_nozero       = find(f(:)   ~= 0 );
fx_nozero      = find(fx(:)  ~= 0 );
fxp_nozero     = find(fxp(:) ~= 0 );
fy_nozero      = find(fy(:)  ~= 0 );
fyp_nozero     = find(fyp(:) ~= 0 );

fid = fopen('fxfyfxpfypf.f90','w');
fprintf(fid,'SUBROUTINE fxfyfxpfypf(Thetta, nf, nfx, nfxp, nfy, nfyp, signout)  \n');

fprintf(fid,['sizef = ', num2str(length(f)),' \n']);
fprintf(fid,['sizex = ', num2str(length(x)),' \n']);
fprintf(fid,['sizey = ', num2str(length(y)),' \n']);

fprintf(fid,['sth  = ', num2str(find( x == 'th_ba1')),' \n']);
fprintf(fid,['seps = ', num2str(find( x == 'eps_ba1')),' \n']);
fprintf(fid,['sgam = ', num2str(find( x == 'gam_ba1')),' \n']);
 
fprintf(fid,['noutput = ', num2str(find( [ y x ] == 'y_cu'     )),' \n']);
fprintf(fid,['npai    = ', num2str(find( [ y x ] == 'pai_cu'   )),' \n']);
fprintf(fid,['nr      = ', num2str(find( [ y x ] == 'r_cu'    )),' \n']);

fprintf( fid,'nf  = 0.0D+0  \n nfx  = 0.0D+00 \n nfxp = 0.0D+00 \n' );
fprintf( fid,'nfy = 0.0D+00 \n nfyp = 0.0D+00 \n' );

for i = 1: length(f_nozero)
    fprintf(fid,['nf(',num2str( f_nozero(i)), ') = ',char(f(f_nozero(i))),'\n']);
end

fprintf(fid,'DO i = 1, sizef \n IF ( ABS(nf(i)) > 1.0D-5) THEN \n');
fprintf(fid,'signout = 1 \n nfx = 0.0D00 \n nfy = 0.0D+00 \n nfxp = 0.0D+00 \n nfyp = 0.0D+00 \n');
fprintf(fid,' write(*,*)  nf \n RETURN \n');
fprintf(fid,' END IF \n END DO \n');

for i = 1: length(fx_nozero)
    fprintf(fid,['nfx(',num2str( fx_nozero(i)), ') = ',char(fx(fx_nozero(i))),' \n ']);
    fprintf(fid,' \n');
end
for i = 1:size(fy_nozero)
    fprintf(fid,['nfy(',num2str( fy_nozero(i)), ') = ',char(fy(fy_nozero(i))),' \n ']);
    fprintf(fid,' \n');
end
for i = 1:size(fxp_nozero)
    fprintf(fid,['nfxp(',num2str( fxp_nozero(i)), ') = ',char(fxp(fxp_nozero(i))),' \n ']);
    fprintf(fid,' \n');
end
for i = 1:size(fyp_nozero)
    fprintf(fid,['nfyp(',num2str( fyp_nozero(i)), ') = ',char(fyp(fyp_nozero(i))),' \n ']);
    fprintf(fid,' \n');
end
fprintf(fid,'END SUBROUTINE fxfyfxpfypf');
fclose(fid);