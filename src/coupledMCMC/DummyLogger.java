package coupledMCMC;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Loggable;


@Description("dummy logger for testing issue #15")
public class DummyLogger extends BEASTObject implements Loggable {

	@Override
	public void initAndValidate() {
	}

	@Override
	public void init(PrintStream out) {
		out.print("dummy\t");
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		out.print("0.0\t");
	}

	@Override
	public void close(PrintStream out) {
	}

}
