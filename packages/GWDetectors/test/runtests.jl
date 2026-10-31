using TestItemRunner

length(ARGS) <= 1 || throw(ArgumentError("expected at most one test selection: default, integration, reference, or all"))
const TEST_SELECTION = isempty(ARGS) ? "default" : only(ARGS)
TEST_SELECTION in ("default", "integration", "reference", "all") ||
    throw(ArgumentError("unknown test selection: $TEST_SELECTION (expected default, integration, reference, or all)"))

function selected_test(item)
    TEST_SELECTION == "all" && return true
    TEST_SELECTION == "default" && return :integration ∉ item.tags
    TEST_SELECTION == "reference" && return :integration in item.tags && :reference in item.tags
    return :integration in item.tags && :reference ∉ item.tags
end

@run_package_tests filter = selected_test
