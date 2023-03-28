
clear
clc

%% path 2 surfaces and electrodes position

path2folder = 'data\head_model_NMTv2\surfaces';

%% plot head model surfaces

for ii=1:4
    
    if ii==1
        disp('Cortex')
        load(fullfile(path2folder, 'tess_cortex_pial.mat'))

        figure;
        V = [1 1 1].*196./255;
        p1 = patch('Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCdata',V, 'FaceAlpha', 0.8); % draw cortex
        set (p1,'EdgeColor','none','FaceColor','flat', ...
            'FaceLighting','phong','diffuseStrength',0.8 );
        camlight right;

    elseif ii==2 % inner skull
        disp('Inner skull')
        load(fullfile(path2folder, 'tess_innerskull_bem.mat'))

        hold on;
        V = [250 153 153]./255;
        p1 = patch('Faces', Faces, 'Vertices', Vertices,'FaceVertexCdata',V, ...
            'FaceAlpha', 0.3); % draw
        set (p1,'EdgeColor','none','FaceColor','flat','FaceLighting', ...
            'phong','diffuseStrength',0.8 );
%         camlight right;
%         lighting phong;
        
    elseif ii==3 % outer skull
        disp('Outer skull')
        load(fullfile(path2folder, 'tess_outerskull_bem.mat'))

        hold on;
        V = [250 204 204]./255;
        p1 = patch('Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCdata',V, 'FaceAlpha', 0.3); % draw
        set (p1,'EdgeColor','none','FaceColor','flat','FaceLighting', ...
            'phong','diffuseStrength',0.8 );
%         camlight right;
%         lighting phong;
        
    else % scalp
        disp('Scalp')
        load(fullfile(path2folder, 'tess_head_bem.mat'))

         hold on;
        V = [250 178 102]./255;
        p1 = patch('Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCdata',V, 'FaceAlpha', 0.2); % draw
        set (p1,'EdgeColor','none','FaceColor','flat','FaceLighting', ...
            'phong','diffuseStrength',0.8 );
        set(gca,'YTickLabel',[],'YTick',[],'YColor',[1 1 1], ...
            'XTickLabel',[],'XTick',[],'XColor',[1 1 1], ...
            'ZTickLabel',[],'ZTick',[],'ZColor',[1 1 1])
%         camlight right;
%         lighting phong;

    end
    
end

%% electrodes

load(fullfile(path2folder,'electrodes.mat'))
selected_electrodes = [6 11 15 31];

hold on;
plot3(electrodes(selected_electrodes, 1), electrodes(selected_electrodes, 2), ...
    electrodes(selected_electrodes, 3), '.y', 'MarkerSize',30)
text(electrodes(selected_electrodes, 1), electrodes(selected_electrodes, 2), ...
    electrodes(selected_electrodes, 3), labels(selected_electrodes))

%% plot dipoles

load(fullfile(path2folder, 'tess_cortex_pial.mat'))
load('SEF_vert_and_orient.mat')

figure;
V = [1 1 1].*196./255; %[0.6350 0.0780 0.1840];
p1 = patch('Faces', Faces, 'Vertices', Vertices,'FaceVertexCdata',V, ...
    'FaceAlpha', 0.8); % draw the outer cortex
set (p1,'EdgeColor','none','FaceColor','flat','FaceLighting', ...
    'phong','diffuseStrength',0.8 );
camlight right;

hold on;
quiver3(Vertices(ind, 1), Vertices(ind, 2), Vertices(ind, 3),...
    orientation(:,1), orientation(:,2), orientation(:,3), 'Color', 'k', ...
    'LineWidth',1.5,'Marker','.', 'MarkerSize',10)

set(gca,'YTickLabel',[],'YTick',[],'YColor',[1 1 1], ...
            'XTickLabel',[],'XTick',[],'XColor',[1 1 1], ...
            'ZTickLabel',[],'ZTick',[],'ZColor',[1 1 1])
