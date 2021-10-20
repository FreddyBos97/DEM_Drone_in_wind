%% Determine the SSE between two states
function SSE = determine_sse(x1,x2,trim)

% Trim of the edges
x1 = x1(trim+1:end-trim);
x2 = x2(trim+1:end-trim);

SSE = sum((x1-x2).^2);
end