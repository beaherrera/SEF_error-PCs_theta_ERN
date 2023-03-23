
%% load fieldtrip defaults settings

ft_defaults

%% initiate brainstorm without the gui

if ~brainstorm('status')
    brainstorm nogui
end