function [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val)
    
    % evaluates the three dominance coefficients for various values of k
    % assuming d = 1

    h1_val = (1/3) / (1/3 + (k_val/(1-k_val)));

    h2_val = 1 / (1 + (k_val/(1-k_val)));

    h3_val = 3 / (3 + (k_val/(1-k_val)));

end