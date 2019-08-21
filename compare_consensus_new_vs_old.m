% load the old one
load 2018-08-01/consensus_neurons.mat
old_neurons = result.consensus_neurons ;
old_neuron_names = result.consensus_neuron_names ;

load 2018-08-01-testing/consensus_neurons.mat
new_neurons = result.consensus_neurons ;
new_neuron_names = result.consensus_neuron_names ;

% ok, there's one more neuron in the new consensus neurons, so that's odd

setdiff(new_neuron_names, old_neuron_names) %  => G-298 is the new one

% ok, looks like Renat added a new neuron since the old runs of this script

all(arrayfun(@(index)(isequal(old_neurons{index}, new_neurons{index})), 1:length(old_neurons)))  % => true

% But that's the only change: the rest of the neurons are the same

% So it has to be a processing issue...
