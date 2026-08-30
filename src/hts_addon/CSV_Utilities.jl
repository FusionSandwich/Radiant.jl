"""Read a simple unquoted comma-separated table used by hash-bound response adapters."""
function _simple_csv_rows(path::AbstractString)
    isfile(path) || error("CSV artifact does not exist: $(path).")
    lines = Tuple{Int,String}[]
    for (line_number,raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line,"#") && continue
        push!(lines,(line_number,line))
    end
    isempty(lines) && error("CSV artifact is empty.")
    header = strip.(split(lines[1][2],',';keepempty=true))
    isempty(header) && error("CSV header is empty.")
    any(isempty,header) && error("CSV headers cannot be empty.")
    length(unique(header)) == length(header) || error("CSV headers must be unique.")
    rows = Dict{String,String}[]
    for (line_number,line) in lines[2:end]
        values = strip.(split(line,',';keepempty=true))
        length(values) == length(header) || error(
            "CSV row $(line_number) has $(length(values)) fields; expected $(length(header)).",
        )
        push!(rows,Dict(header[index] => values[index] for index in eachindex(header)))
    end
    isempty(rows) && error("CSV artifact contains no data rows.")
    return header,rows
end

function _csv_float(row::AbstractDict,key::AbstractString)
    key_string = String(key)
    haskey(row,key_string) || error("CSV row is missing column $(key_string).")
    text = strip(string(row[key_string]))
    isempty(text) && error("CSV value for $(key_string) is empty.")
    value = try
        parse(Float64,text)
    catch
        error("CSV value for $(key_string) is not numeric: $(text).")
    end
    isfinite(value) || error("CSV value for $(key_string) must be finite.")
    return value
end

function _csv_int(row::AbstractDict,key::AbstractString)
    value = _csv_float(row,key)
    isinteger(value) || error("CSV value for $(key) must be an integer.")
    return Int64(round(value))
end

function _csv_required_columns(header::AbstractVector,required::AbstractVector)
    missing = [String(key) for key in required if String(key) ∉ header]
    isempty(missing) || error("CSV artifact is missing columns: $(join(missing,',' )).")
    return true
end
