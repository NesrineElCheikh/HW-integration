module HW_int

			using FastGaussQuadrature
			using Roots
			using Sobol
			using PyPlot
			using Distributions
			using DataFrames

			function question_1b(n)
				3*sum(gausslegendre(n)[2][i]*(1.5*gausslegendre(n)[1][i]+2.5)^(-0.5) for i=1:n)
			end

			function question_1c(n)
				A = Float64[]
					for i=1:n
					push!(A, (1:4))
					end
				3*(1/n)*sum(2*(A[i])^(-0.5) for i=1:n)
			end

			function question_1d(n)
				s = SobolSeq(1)
				B = []
					for i=1:n
					push!(B, 1+(next(s))*3)
					end
				3*(2/n)*sum((B[i][1])^(-0.5) for i=1:n)
			end

			function question_2a(n)
				sigma = [0.01 0.01; 0.01 0.02]
				OT = chol(sigma)
				O = transpose(OT)
				C = Any[]
					push!(C,repeat(gausshermite(n)[1],inner=[1],outer=[n]))
					push!(C,repeat(gausshermite(n)[1],inner=[n],outer=[1]))
				weights =(2*pi)^(-1)*det(sigma)^(-0.5)*kron(gausshermite(n)[2],gausshermite(n)[2])
				prenodes = [O[1,1]*C[1, :] + O[1,2]*C[2, :]; O[2,1]*C[1, :]+O[2,2]*C[2, :]]
				nodes = exp.(prenodes)
				df = hcat(weights, nodes[1, :][1], nodes[2, :][1])
				X(p) = sum(df[:, 1] .* (df[:, 2]/p .+ df[:, 3]/sqrt(p) - 2))
				exp_p = fzero(X, 0, Inf)
				gh_est = [exp_p]
			end


					# function to run all questions
			function runall(n)
				println("running all questions of HW-integration:")
				println("results of question 1")
				q1b = question_1b(n)	# make sure your function prints some kind of result!
				q1c = question_1c(n)
				q1d = question_1d(n)
				println("Question 1b")
				println(q1b)
				println("Question 1c")
				println(q1c)
				println("Question 1d")
				println(q1d)
				println("")
				println("results of question 2 (expectation, variance):")
				q2 = question_2a(n)
				println("Question 2a")
				println(q2)
				q2b = question_2b(n)
				println("Question 2b")
				println(q2b)
				println("end of HW-integration")
			end
		end
