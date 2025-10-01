function plot_mode_Efields(X,Y,E_vec,cmap,h_params,field_true)

Ex=E_vec(:,:,1);
Ey=E_vec(:,:,2);
Ez=E_vec(:,:,3);

norm_fac_x=max(max(abs(Ex).^2));
norm_fac_y=max(max(abs(Ey).^2));
norm_fac_z=max(max(abs(Ez).^2));
norm_fac_e=max([norm_fac_x,norm_fac_y,norm_fac_z]);

h1=h_params.h1;
h2=h_params.h2;
h3=h_params.h3;
if isfield(h_params,'w_r')
    w_r=h_params.w_r;
    w_b=h_params.w_b;
    gap=h_params.gap;
    x_crop_start=-(w_b/2+gap+w_r+500);
    x_crop_stop=w_b+1000;
%     x_crop_start=-w_r/2-1000;
%     x_crop_stop=w_r/2+gap+w_b+500;
    x_crop_len=x_crop_stop-x_crop_start;
    
    y_cent=(h1+h2+h3)-h1-0.5*h2;
    y_crop_top=y_cent+0.5*x_crop_len;
    y_crop_bottom=y_cent-0.5*x_crop_len;
    
else
    w=h_params.w;
    y_crop_top=(h1+h2+h3)-h1+1000;
    y_crop_bottom=y_crop_top-h2-2000;
    x_crop_fac=abs(y_crop_top-y_crop_bottom)/w/2;
    x_crop_start=-1*x_crop_fac*w;
    x_crop_stop=1*x_crop_fac*w;
    
end

if ~field_true
    figure;
    subplot(1,3,1);
    pcolor(X'*10^-3,Y'*10^-3,abs(Ex).^2/norm_fac_e); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('|E_x|^2'); set(gca,'fontsize',14);
    subplot(1,3,2);
    pcolor(X'*10^-3,Y'*10^-3,abs(Ey).^2/norm_fac_e); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('|E_y|^2'); set(gca,'fontsize',14);
    xlabel('x (\mum)');
    subplot(1,3,3);
    pcolor(X'*10^-3,Y'*10^-3,abs(Ez).^2/norm_fac_e); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('|E_z|^2'); set(gca,'fontsize',14);
    xlabel('x (\mum)');
else
    figure;
    subplot(1,3,1);
    pcolor(X'*10^-3,Y'*10^-3,real(Ex)); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('Re(E_x)'); set(gca,'fontsize',14);
    subplot(1,3,2);
    pcolor(X'*10^-3,Y'*10^-3,real(Ey)); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('Re(E_y)'); set(gca,'fontsize',14);
    xlabel('x (\mum)');
    subplot(1,3,3);
    pcolor(X'*10^-3,Y'*10^-3,real(Ez)); view(2); shading interp
    xlim([x_crop_start*10^-3 x_crop_stop*10^-3]); ylim([y_crop_bottom*10^-3 y_crop_top*10^-3]);
    axis square; caxis manual; caxis([0 1]); colormap(cmap)
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title('Re(E_z)'); set(gca,'fontsize',14);
    xlabel('x (\mum)');
end