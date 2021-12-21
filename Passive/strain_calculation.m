function [strain_images] = strain_calculation(displacement,pixel_size,refrac,ax)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

pixel_res=pixel_size/refrac;
%ax=50;
z_kern=ax;

strain_images=zeros(size(displacement));
%strain_images(:,:,((tf-t0)+1):end)=[];

for n=1:size(displacement,3)
    
    displacement_frame=medfilt2(displacement(:,:,n)',[7,7]);
    % displacement=displacement';
    [bscanwidth,how_many_depths] = size(displacement_frame);
    b_scan_strain = single(zeros(bscanwidth,how_many_depths-z_kern));
    
    for z1 = 1:how_many_depths-z_kern+1
        sub_displacement = squeeze(displacement_frame(:,z1:z1+z_kern-1))';
        % Vandermonde matrix
        sub_strain = ( [(1:size(sub_displacement))',ones(size(sub_displacement,1),1) ] \ sub_displacement)./pixel_res;
%         b_scan_strain(:,z1:z1+z_kern-1) = 1000.*repmat(sub_strain(1,:)',1,z_kern);
                b_scan_strain(:,z1:z1+z_kern-1) = repmat(sub_strain(1,:)',1,z_kern);

        %         figure(4);
        %         set(gcf,'units','normalized','outerposition',[0.5 0 0.5 0.5]);
        %         imagesc(strain_cpu');
        %         cbar1 = colorbar; colormap(jet); caxis([0,5]);
        %         ylabel(cbar1,'Strain (m\epsilon)','Rotation',270,'Position',[5,2.5]);
        %         set(gca,'XTick',x_ticks);
        %         set(gca,'XTickLabels',x_tick_labels);
        %         xlabel('X (mm)');
        %         set(gca,'YTick',y_ticks);
        %         set(gca,'YTickLabels',y_tick_labels);
        %         ylabel('Depth (mm)');
        %         set(gca,'FontSize',20);
        %         title({'Strain';pathname},'Interpreter','none','FontSize',12);
        %
        %         wait_bar = waitbar(z1/(how_many_depths-z_kern),wait_bar,...
        %             ['Calculating Strain (',num2str(z1/(how_many_depths-z_kern)*100),'%)']);
    end
    %m=(n)+1;
    
    strain_images(:,:,n)=b_scan_strain';
end
end

