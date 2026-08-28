# To use Python-like notation for methods.
const RadiantObject = Union{
    Material,
    Cross_Sections,
    Multigroup_Cross_Sections,
    Geometry,
    SN,
    Solvers,
    Surface_Source,
    Volume_Source,
    Fixed_Sources,
    Source,
    Sources,
    Computation_Unit,
    Flux_Per_Particle,
    Flux,
    Interaction,
    Particle,
    GN,
    CP,
    Electromagnetic_Field,
    Source_Normalization,
    Boundary_Angular_Current_Source,
    Anisotropic_Volume_Source,
    Transport_Balance,
    Transport_Ownership_Record,
    Transport_Ownership_Map,
    Layer_Definition,
    Tape_Stack_Definition,
    Energy_Partition,
    Fibre_Radial_Layer,
    Fibre_Diagnostic_Definition,
    Detector_Converter_Layer,
    Detector_Thermal_Model,
    HTS_Detector_Definition,
    Physics_Coverage_Record,
    Physics_Coverage_Register,
    Protected_Response,
    Protected_Response_Registry
}

"""
    Base.getproperty(object::RadiantObject, s::Symbol)

Enable calling a method `m(this::T,...)` using `object.m(...)` while preserving ordinary field
access. This compatibility layer is retained for the existing Radiant user API.
"""
function Base.getproperty(object::RadiantObject,s::Symbol)
    if hasfield(typeof(object),s)
        return getfield(object,s)
    end
    if isdefined(Radiant,s) && getproperty(Radiant,s) isa Function
        return function(arguments...)
            function_object = getproperty(Radiant,s)
            return function_object(object,arguments...)
        end
    end
    type = typeof(object)
    error("The object of type $(type) has no field or Radiant method $(s).")
end
