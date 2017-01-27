process = qp2qpg qP2qPg qQ2Ppg qg2qpP    qq2qqg qQ2qQg qg2qqQ    qQ2ggg qg2qgg Qg2Qgg gg2qQg    gg2ggg

all: $(process)

$(process):
	form -q -D $@ me2to3.frm
	@mv me2to3.temp.m $@.m

clean:
	rm *.m
