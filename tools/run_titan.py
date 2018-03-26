from wraprun import Wraprun


def main():

    systems = ['e255k', 'e255v']
    num_replicas = 25

    job = Wraprun()

    for system in systems:
        for replica in range(num_replicas):
            job.add_task("-n 1 serial $MINICONDA3_PYTHON run_esmacs.py --system {} --replica {}".format(system, num_replicas))

    job.launch()


if __name__ == '__main__':
    main()
