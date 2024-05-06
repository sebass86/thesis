using Distributed

function test_processors(p_count)
    rmprocs(procs())
    addprocs(p_count)
    @sync @distributed for i in 1:100_000  # Assuming each iteration is a stand-alone task
        # Your computation here
    end
end

# Run tests
for p in [2, 4, 8, 12, 30]  # Adjust based on your hardware
    println("Testing with $p processors.")
    @time test_processors(p)
end
rmprocs(workers())