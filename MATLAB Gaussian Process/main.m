%Running this script will add the necessary files to the path and open a
%GUI with the tests done 
%The SONIG version menu refers to how the SONIG was trained

%Trained on different class - was trained on an entirely different set of
%trajectories this one gives numerical problems whenever it is used on the
%the different class

%Really similiar - Some handpicked trajectories that are pretty close to
%each other

%Same class, all trajectories - trained on almost every trajectory of the
%same class, except for a few that are used for validation and test

%Same class, no time as a feature - wanted to show what happens when I
%remove time as a feature, doesn't really matter what I do it always looks
%worse than when I include it. Was trained on most trajectories, same as
%above

%Similiar trajectories, optimized - this one was trained on the similiar
%trajectories chosen mostly because it is a smaller set so it doesn't take
%forever to find the parameters


%I also tried to train without time as feature with the different class,
%resulting in staggeringly ugly results that no innocent human should be ever
%subject to
%Just kidding but it didn't look any good

%The test trajectory menu refers to which trajectory of the test set is
%used

%The number of steps predicted can be changed on the editable box below

%The current position can be moved by pressing the buttons

addpath(genpath('SONIG'))
addpath(genpath('scripts'))

gui

