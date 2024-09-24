# Project     : merge_fastq
# File Name   : bsub.py
# Description : Methods for submitting LSF jobs at WashU.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Sep 24 15:33:28 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import os
import yaml  # type: ignore
from datetime import datetime
from shutil import copyfile


class Bsub:
    """
    A class for executing bsub jobs at WashU.

    This class is used for creating, tracking, and executing LSF bsub
    jobs on the WashU compute1 compute cluster. This code is most useful
    when writing scripts to loop through a number of independent tasks
    that should receive their own bsub job.

    Imports
    -------

    os
    yaml
    datetime
    shutil (copyfile)

    Parameters
    ----------

    docker_volumes : dict
        Dictionary holds key:value pairs for Docker-style volume
        mappings.

    docker_preserve_environment : bool, default: False
        Turns preserving the local RIS environment on/off when running a
        Docker instance.

    docker_image : str
        Sets the Docker image from which to pull --- e.g. DockerHub
        image.

    memory_max : str, default: 8G
        Sets maximum memory for a bsub job.

    group : str
        Sets the computation group for a bsub job.

    queue : str
        Sets the computation queue for a bsub job.

    error_log : str, default: bsub.err
        Sets the file name for writing a bsub job's error messages
        (STDERR).

    output_log : str, default: bsub.out
        Sets the file name for writing a bsub job's output (STDOUT).

    config : str, default: config.yaml
        Sets the file name of the bsub job's YAML configuration file.

    bsub_command_name : str, default: bsub_cmd.sh
        Sets the file name of the bsub command.

    command_name : str, default: cmd.sh
        Sets the file name of the shell command run by the bsub command.

    log_dir : str, default: __bsub
        Sets the name of the bsub log directory, where log files will be
        written.

    job_name : str
        Set the name of the bsub job.

    kill_time : str
        Sets the runtime limit (minutes) of the job.

    number_of_tasks : str
        Submits a parallel job and specifies the number of tasks in the
        job.

    email : str
        Email address to notify user of bsub job statuses.

    resource_memory : str, default: 8G
        Sets the bsub job resource memory requirement.

    resource_tmp : str
        Sets bsub job's temporary memory requirement.

    resource_usage_memory : str, default: 8G
        Sets the bsub job's resource usage memory.

    resource_usage_tmp : str
        Sets the bsub job's resource usage temporary memory.

    resource_span_hosts : int, default: 1
        Sets the number of hosts to request for the bsub job.

    current_working_dir : str, default: os.getcwd()
        Sets the current working directory for the bsub jobs. Unless
        specified, the default is the current working directory where
        the script calling the class was executed.

    command : str | list
        If command is set to a string (command file path) then the code
        will try to ruin the provided file. If command is set to a list,
        the list is taken as a literal list of commands and a new
        command file will be written based on the list of commands.

    Methods
    -------

    execute(self, dry: bool = False) -> None
        For a given bsub job object, executes the instance and runs the
        job on the LSF system. If dry=True, then all bsub related files
        will be written to disk, but the job will not be executed on the
        LSF system; useful for assessment prior to running a job.

    print_self(self) -> None
        Dumps the contents of the object's instance in YAML format to
        STDOUT.

    print_bsub_command(self) -> None
        Prints the object's full bsub command to STDOUT.

    reset_execution_counter(self) -> None
        The class keeps track of the ids of all bsub jobs (objects)
        created. If you want to have multiple output sessions, with
        their own ids, under the same class import, the class execution
        counter will need to be reset for each session by calling this
        method prior to a new session.

    The following methods set object variables. See Attributes above for
    variable descriptions.

    set_bsub_command_name(self, value: str) -> None
    set_command(self, value: str) -> None
    set_command_name(self, value: str) -> None
    set_config(self, value: str) -> None
    set_current_working_directory(self, value: str) -> None
    set_docker_image(self, value: str) -> None
    set_docker_preserve_environment(self, value: bool) -> None
    set_docker_volumes(self, value: dict) -> None
    set_email(self, value: str) -> None
    set_error_log(self, value: str) -> None
    set_group(self, value: str) -> None
    set_job_name(self, value: str) -> None
    set_kill_time(self, value: str) -> None
    set_log_dir(self, value: str) -> None
    set_memory_max(self, value: str) -> None
    set_number_of_tasks(self, value: str) -> None
    set_output_log(self, value: str) -> None
    set_queue(self, value: str) -> None
    set_resource_memory(self, value: str) -> None
    set_resource_span_hosts(self, value: int) -> None
    set_resource_tmp(self, value: str) -> None
    set_resource_usage_memory(self, value: str) -> None
    set_resource_usage_tmp(self, value: str) -> None

    Raises
    ------

    TypeError
        ERROR MESSAGE: The setter encountered a disparate data type
        than the one expected.

    ValueError
        ERROR MESSAGE: A required value was not found.

    FileExistsError
        The log_dir directory already existed.

    Examples
    --------

    A very simple example of running a single bsub job by passing a list
    of shell commands to a bsub object. Note that the job writes all of
    the bsub files but does not actually run the job on LSF --- i.e.
    dry=True.

    >>> import sys
    >>> sys.path.append('.')
    >>> from washu.ris.bsub import Bsub
    >>>
    >>> cmd = ['print("Hello World!")', '# __END__']
    >>>
    >>> job = Bsub(
    >>>     docker_image='twylie/bfx_toolbox',
    >>>     group='twylie_group',
    >>>     queue='research',
    >>>     command=cmd
    >>> )
    >>>
    >>> job.execute(dry=True)

    A more complex example showing creating multiple bsub objects within
    a loop. Here, we have set more object variables.

    >>> from washu.ris.bsub import Bsub
    >>>
    >>> volumes = (
    >>> {'/storage1/fs1/twylie/Active': '/storage1/fs1/twylie/Active'}
    >>> )
    >>>
    >>> for i in range(0, 10):
    >>>
    >>>     cmd = ['echo "Hello World number {}!"'.format(i), '# __END__']
    >>>
    >>>     error_log = '{}_bsub.err'.format(i)
    >>>     command_name = '{}_cmd.sh'.format(i)
    >>>     job_name = 'job_{}'.format(i)
    >>>     output_log = '{}_bsub.out'.format(i)
    >>>     config = '{}_config.yaml'.format(i)
    >>>     bsub_command_name = '{}_bsub_cmd.sh'.format(i)
    >>>
    >>>     job = Bsub(
    >>>         docker_volumes=volumes,
    >>>         docker_image='twylie/bfx_toolbox',
    >>>         group='twylie_group',
    >>>         queue='research',
    >>>         command=cmd,
    >>>         error_log=error_log,
    >>>         output_log=output_log,
    >>>         command_name=command_name,
    >>>         config=config,
    >>>         bsub_command_name=bsub_command_name,
    >>>         number_of_tasks='1',
    >>>         kill_time='120',
    >>>         email='twylie@wustl.edu',
    >>>         job_name=job_name,
    >>>         resource_tmp='2G',
    >>>         resource_span_hosts=10,
    >>>         resource_usage_memory='20G',
    >>>         resource_usage_tmp='2G'
    >>>     )
    >>>
    >>>     job.execute()
    """

    execution_counter = 0  # Counts execute() instances.

    def __init__(
            self,
            docker_volumes: dict = {},
            docker_preserve_environment: bool = False,
            docker_image: str = '',
            memory_max: str = '8G',
            group: str = '',
            queue: str = '',
            error_log: str = 'bsub.err',
            output_log: str = 'bsub.out',
            config: str = 'config.yaml',
            bsub_command_name: str = 'bsub_cmd.sh',
            command_name: str = 'cmd.sh',
            log_dir: str = '__bsub',
            job_name: str = '',
            kill_time: str = '',
            number_of_tasks: str = '',
            email: str = '',
            resource_memory: str = '8G',
            resource_tmp: str = '',
            resource_usage_memory: str = '8G',
            resource_usage_tmp: str = '',
            resource_span_hosts: int = 1,
            current_working_dir: str = os.getcwd(),
            command: str = ''
    ) -> None:
        """
        Construct the class.

        The parameters collected here are outlined above under the Bsub
        class Parameters section.
        """

        self.set_docker_volumes(docker_volumes)
        self.set_docker_preserve_environment(docker_preserve_environment)
        self.set_docker_image(docker_image)
        self.set_memory_max(memory_max)
        self.set_group(group)
        self.set_queue(queue)
        self.set_error_log(error_log)
        self.set_output_log(output_log)
        self.set_config(config)
        self.set_bsub_command_name(bsub_command_name)
        self.set_command_name(command_name)
        self.set_log_dir(log_dir)
        self.set_job_name(job_name)
        self.set_kill_time(kill_time)
        self.set_number_of_tasks(number_of_tasks)
        self.set_email(email)
        self.set_resource_memory(resource_memory)
        self.set_resource_tmp(resource_tmp)
        self.set_resource_usage_memory(resource_usage_memory)
        self.set_resource_usage_tmp(resource_usage_tmp)
        self.set_resource_span_hosts(resource_span_hosts)
        self.set_current_working_directory(current_working_dir)
        self.set_command(command)

        return

    # SETTER METHODS ##########################################################

    def set_docker_volumes(self, value: dict) -> None:
        if type(value) is not dict:
            self.__setter_error(
                method='set_docker_volumes',
                error_code='must be of dict type'
            )
        self.docker_volumes = value
        return

    def set_docker_preserve_environment(self, value: bool) -> None:
        if type(value) is not bool:
            self.__setter_error(
                method='set_docker_preserve_environment',
                error_code='must be of bool type'
            )
        self.docker_preserve_environment = value
        return

    def set_memory_max(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_memory_limit',
                error_code='must be of str type'
            )
        self.memory_max = value
        return

    def set_docker_image(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_docker_image',
                error_code='must be of str type'
            )
        self.docker_image = value
        return

    def set_group(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_group',
                error_code='must be of str type'
            )
        self.group = value
        return

    def set_queue(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_queue',
                error_code='must be of str type'
            )
        self.queue = value
        return

    def set_error_log(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_error_log',
                error_code='must be of str type'
            )
        self.error_log = value
        return

    def set_output_log(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_output_log',
                error_code='must be of str type'
            )
        self.output_log = value
        return

    def set_config(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_config',
                error_code='must be of str type'
            )
        self.config = value
        return

    def set_bsub_command_name(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_bsub_command_name',
                error_code='must be of str type'
            )
        self.bsub_command_name = value
        return

    def set_command_name(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_command_name',
                error_code='must be of str type'
            )
        self.command_name = value
        return

    def set_log_dir(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_log_dir',
                error_code='must be of str type'
            )
        self.log_dir = value
        return

    def set_job_name(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_job_name',
                error_code='must be of str type'
            )
        self.job_name = value
        return

    def set_kill_time(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_kill_time',
                error_code='must be of str type'
            )
        self.kill_time = value
        return

    def set_number_of_tasks(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_number_of_tasks',
                error_code='must be of str type'
            )
        self.number_of_tasks = value
        return

    def set_email(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_email',
                error_code='must be of str type'
            )
        self.email = value
        return

    def set_resource_memory(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_resource_memory',
                error_code='must be of str type'
            )
        self.resource_memory = value
        return

    def set_resource_tmp(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_resource_tmp',
                error_code='must be of str type'
            )
        self.resource_tmp = value
        return

    def set_resource_usage_memory(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_resource_usage_memory',
                error_code='must be of str type'
            )
        self.resource_usage_memory = value
        return

    def set_resource_usage_tmp(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_resource_usage_tmp',
                error_code='must be of str type'
            )
        self.resource_usage_tmp = value
        return

    def set_resource_span_hosts(self, value: int) -> None:
        if type(value) is not int:
            self.__setter_error(
                method='set_resource_span_hosts',
                error_code='must be of int type'
            )
        self.resource_span_hosts = value
        return

    def set_current_working_directory(self, value: str) -> None:
        if type(value) is not str:
            self.__setter_error(
                method='set_current_working_directory',
                error_code='must be of str type'
            )
        self.current_working_dir = value
        return

    def set_command(self, value: str) -> None:
        if type(value) is not str and type(value) is not list:
            self.__setter_error(
                method='set_command',
                error_code='must be of str or list type'
            )
        self.command = value
        return

    # PRIVATE: ERROR MESSAGES #################################################

    def __setter_error(self, method: str, error_code: str) -> None:
        """
        Handle a setter error.

        A generic way of handling an error encountered when trying to
        set on of the object's attributes.

        Parameters
        ----------
        method : str
            States the method from where the error orginated.
        error_code : str
            The error code.

        Raises
        ------
        TypeError
            The setter encountered a disparate data type than the one
            expected.

        Examples
        --------
        self.__setter_error(
            method='set_command',
            error_code='must be of str or list type'
        )
        """
        msg = f'method: {method} error: {error_code}'
        raise TypeError(msg)

    def __required_error(self, variable: str, error_code: str) -> None:
        """
        Handle a required attribute not being set.

        A generic way of handling an error message when a required
        attribute is missing from the object.

        Parameters
        ----------
        variable : str
            The variable being evaluated.
        error_code : str
            The variable being evaluated.

        Raises
        ------
        ValueError
            A required value was not found.

        Examples
        --------
        self.__required_error(
            variable='docker_image',
            error_code='is empty, but must be defined'
        )
        """
        msg = f'variable: {variable} error: {error_code}'
        raise ValueError(msg)

    # PRIVATE: EVAL SETTINGS ##################################################

    def __eval_settings(self) -> None:
        """
        Check object for required attributes.

        REQUIRED
        --------
        docker_image
        group
        queue
        command

        Examples
        --------
        self.__eval_settings()
        """

        if not self.docker_image:
            self.__required_error(
                variable='docker_image',
                error_code='is empty, but must be defined'
            )

        if not self.group:
            self.__required_error(
                variable='group',
                error_code='is empty, but must be defined'
            )

        if not self.queue:
            self.__required_error(
                variable='queue',
                error_code='is empty, but must be defined'
            )

        if not self.command:
            self.__required_error(
                variable='command',
                error_code='is empty, but must be defined'
            )

        return

    # PRIVATE: FORMULATE COMMAND ##############################################

    def __formulate_bsub_command(self) -> None:
        """
        Create the object's bsub commands.

        Examples
        --------
        self.__formulate_bsub_command()
        """

        self.__eval_settings()

        log_dir = os.path.join(self.current_working_dir, self.log_dir)

        cmd_error_log = '-e {}'.format(os.path.join(log_dir, self.error_log))

        cmd_output_log = '-o {}'.format(os.path.join(log_dir, self.output_log))

        self.docker_volumes.update(
            {self.current_working_dir: self.current_working_dir}
        )

        volume_pairs = list()
        for i in self.docker_volumes.keys():
            pair = ':'.join([i, self.docker_volumes[i]])
            volume_pairs.append(pair)

        cmd_docker_volumes = 'LSF_DOCKER_VOLUMES="{}"'.format(
            ' '.join(volume_pairs)
        )

        cmd_docker_preserve_environment = (
            'LSF_DOCKER_PRESERVE_ENVIRONMENT={}'.format(
                str(self.docker_preserve_environment).lower()
            )
        )

        cmd_group = '-G {}'.format(self.group)

        cmd_queue = '-q {}'.format(self.queue)

        cmd_memory_max = '-M {}'.format(self.memory_max)

        cmd_docker_image = "-a 'docker({})'".format(self.docker_image)

        if self.resource_tmp:
            resource_memory = 'select[mem>{} && tmp>{}]'.format(
                self.resource_memory, self.resource_tmp
            )
        else:
            resource_memory = 'select[mem>{}]'.format(self.resource_memory)

        resource_span_hosts = 'span[hosts={}]'.format(self.resource_span_hosts)

        if self.resource_usage_tmp:
            resource_usage_memory = 'rusage[mem={}, tmp={}]'.format(
                self.resource_usage_memory, self.resource_usage_tmp
            )
        else:
            resource_usage_memory = (
                'rusage[mem={}]'.format(self.resource_usage_memory)
            )

        resource_params = ' '.join([
            resource_memory,
            resource_span_hosts,
            resource_usage_memory
        ])

        resource_command = '-R "{}"'.format(resource_params)

        command_file = os.path.join(
            self.current_working_dir, self.log_dir, self.command_name
        )
        self.command_file = command_file
        cmd = 'sh {}'.format(self.command_file)

        bsub_cmd_head = ' '.join([
            cmd_docker_preserve_environment,
            cmd_docker_volumes,
            'bsub',
        ])

        bsub_cmd_tail = ' '.join([
            resource_command,
            cmd_memory_max,
            cmd_group,
            cmd_queue,
            cmd_output_log,
            cmd_error_log,
            cmd_docker_image,
            cmd
        ])

        optional_params = list()

        if self.job_name:
            optional_job_name = '-J "{}"'.format(self.job_name)
            optional_params.append(optional_job_name)

        if self.number_of_tasks:
            optinal_number_of_tasks = '-n {}'.format(self.number_of_tasks)
            optional_params.append(optinal_number_of_tasks)

        if self.kill_time:
            optional_kill_time = '-W {}'.format(self.kill_time)
            optional_params.append(optional_kill_time)

        if self.email:
            email = '-N -u "{}"'.format(self.email)
            optional_params.append(email)

        if optional_params:
            j_optional_params = ' '.join(optional_params)
            self.full_bsub_command = ' '.join(
                [bsub_cmd_head, j_optional_params, bsub_cmd_tail]
            )
        else:
            self.full_bsub_command = ' '.join([bsub_cmd_head, bsub_cmd_tail])

        now = datetime.now()
        self.date_time = now.strftime('%d/%m/%Y %I:%M:%S %p')

        return

    # PRIVATE: WRITE CONFIG FILE ##############################################

    def __write_config_file(self) -> None:
        """
        Write the config file to disk.

        Examples
        --------
        self.__write_config_file()
        """
        config_file = os.path.join(self.log_dir, self.config)
        with open(config_file, 'w') as fh:
            yaml.dump(self.__dict__, fh)
        return

    # PRIVATE: MAKE LOG DIR ###################################################

    def __make_log_dir(self) -> None:
        """
        Make the log directory.

        Raises
        ------
        FileExistsError
            The log_dir directory already existed.

        Examples
        --------
        self.__make_log_dir()
        """

        log_dir = os.path.join(self.current_working_dir, self.log_dir)

        if os.path.isdir(log_dir):
            raise FileExistsError(
                'The log_dir directory already exists.', log_dir
            )
        else:
            os.mkdir(log_dir)

        return

    # PRIVATE: WRITE BSUB COMMAND FILE ########################################

    def __write_bsub_command_file(self) -> None:
        """
        Write the bsub command file to disk.

        Examples
        --------
        self.__write_bsub_command_file()
        """

        bsub_command_path = os.path.join(
            self.current_working_dir, self.log_dir, self.bsub_command_name
        )
        self.bsub_command_file = bsub_command_path

        with open(bsub_command_path, 'w') as fh:
            fh.write(self.full_bsub_command + '\n')

        return

    # PRIVATE: WRITE COMMAND FILE #############################################

    def __write_command_file(self) -> None:
        """
        Write the shell command file to disk.

        Examples
        --------
        self.__write_command_file()
        """

        if type(self.command) is str:
            # Command file already exists.
            copyfile(self.command, self.command_file)
        elif type(self.command) is list:
            # We will write the command file.
            with open(self.command_file, 'w') as fh:
                for line in self.command:
                    fh.write(line + '\n')

        return

    # EXECUTE #################################################################

    def execute(self, dry: bool = False) -> None:
        """
        Execute the bsub object's job.

        This is the main method for executing a bsub object's job. There
        is an option for dry run execution or actually submitting the
        bsub job to the LSF system.

        Parameters
        ----------
        dry : bool, default: False
            Either dry run the job or submit it to the LSF system.

        Examples
        --------
        job.execute(dry=True)
        """

        self.__formulate_bsub_command()

        if os.path.isdir(self.log_dir) is False:
            self.__make_log_dir()

        self.__write_config_file()

        self.__write_bsub_command_file()

        self.__write_command_file()

        Bsub.execution_counter += 1

        if dry is True:
            pass
        else:
            print('Running job {}: {}'.format(
                Bsub.execution_counter, self.bsub_command_file)
                  )
            os.system('sh ' + self.bsub_command_file)

        return

    # PRINT CONFIG ############################################################

    def print_self(self) -> None:
        """
        Print a YAML version of the config to STDOUT.

        Examples
        --------
        job.print_self()
        """
        self.__formulate_bsub_command()
        print('')
        print(yaml.dump(self, sort_keys=False))
        return

    # PRINT BSUB COMMAND ######################################################

    def print_bsub_command(self) -> None:
        """
        Print bsub command to STDOUT.

        Examples
        --------
        job.print_bsub_command()
        """
        self.__formulate_bsub_command()
        print(self.full_bsub_command)
        return

    # RESET EXECUTION COUNTER #################################################

    def reset_execution_counter(self) -> None:
        """
        Reset the class variable execution counter.

        The class keeps track of the ids of all bsub jobs (objects)
        created. If you want to have multiple output sessions, with
        their own ids, under the same class import, the class execution
        counter will need to be reset for each session by calling this
        method prior to a new session.

        Examples
        --------
        job.reset_execution_counter()
        """
        Bsub.execution_counter = 0
        return

# __END__
